import numpy as np
import random
import uuid


class Pool:

    # this could actually just be a dict, pointing to the member
    # cells. the UID is not really necessary, except for displaying
    # and then it's also only necessary to display something
    # to distinguish pools from one another.
    # otherwise, this could just be a dict with a few values.
    # which, if it were that, this whole thing would be more viable
    # for multiprocessing.

    def __init__(self, waterlevelheight):
        self.inflow = 0
        self.outflow = 0
        self.averaged_total_height = 0
        self.outflow_points = []
        self.member_cells = []
        self.uid = str(uuid.uuid4())
        # if some river points ot the outside, drain a lot of water
        # from all pool cells, to lower the level.
        self.riverdrain = False
        self.equalize_elevation = True

    def equalize_main(self, meta_info):

        total_rml = []
        rml = []

        while True:
            if rml != []:
                total_rml += rml
            rml = []

            if len(self.member_cells) == 0:
                return []

            self.remove_empty_get_average(meta_info, rml)
            self.remove_outliers(meta_info, rml)

            if rml != []:
                for x in rml:
                    total_rml.append(x)
                    self.member_cells.remove(x)
                continue

            self.final_equalize(meta_info)

            break
        for x in rml:
            self.member_cells.remove(x)

        return total_rml

    def remove_empty_get_average(self, meta_info, rml):
        total_el = 0

        for cell in self.member_cells:
            wl = meta_info[cell]["waterlevel"]
            el = meta_info[cell]["elevation"]

            # if there is no water in this one, it's definitely not part of the pool.
            # don't even count it here.
            if wl == 0:
                rml.append(cell)
                continue
            total_el += (wl + el)

        self.averaged_total_height = total_el / len(self.member_cells)

    def remove_outliers(self, meta_info, rml):
        # figure out if it's realistic to see this one as part of my pool.
        for cell in self.member_cells:

            wl = meta_info[cell]["waterlevel"]
            el = meta_info[cell]["elevation"]

            if cell in rml:
                continue

            my_sum = (wl+el)
            my_diff = self.averaged_total_height - my_sum

            # if it's flowing out
            if my_diff < 0 and abs(my_diff) > meta_info[cell]["waterlevel"]:
                rml.append(cell)
                continue

    def equalize_elevation_member_cells(self, meta_info):
        el_sum = 0
        for cell in self.member_cells:
            el = meta_info[cell]["elevation"]
            el_sum += el
        memberl = len(self.member_cells)

        other_avg = (el_sum/memberl) * 0.85

        for cell in self.member_cells:
            this_cell = meta_info[cell]["elevation"] * 0.15
            meta_info[cell]["elevation"] = other_avg + this_cell

    def final_equalize(self, meta_info):

        # if I want a good drain map, I may have to equalize this.
        equalize_elevation = True
        if self.equalize_elevation:
            self.equalize_elevation_member_cells(meta_info)

        for cell in self.member_cells:
            wl = meta_info[cell]["waterlevel"]
            el = meta_info[cell]["elevation"]
            my_sum = (wl+el)
            my_diff = self.averaged_total_height - my_sum

            if self.riverdrain:
                # meta_info[cell]["waterlevel"] *= 0.65
                # meta_info[cell]["elevation"] *= 0.9
                self.riverdrain = False

            meta_info[cell]["waterlevel"] += my_diff
            meta_info[cell]["waterlevel"] = max(
                meta_info[cell]["waterlevel"], 0)


def generate_basic_data():

    xlen, ylen = 15, 15

    centers = []
    x = 0.5
    while x < xlen:
        y = 0.5
        while y < ylen:
            centers.append((x, y, 0))
            y += 1
        x += 1

    all_neighbors = {}
    total_counter = 0
    x = 0.5
    while x < xlen:
        y = 0.5
        while y < ylen:
            # see if these exist
            up = (x, y+1, 0)
            down = (x, y-1, 0)
            left = (x-1, y, 0)
            right = (x+1, y, 0)
            directions = [up, down, left, right]
            neighbors = []
            for p in directions:
                if p in centers:
                    my_index = centers.index(p)
                    neighbors.append(my_index)
            all_neighbors[total_counter] = neighbors
            total_counter += 1
            y += 1
        x += 1

    meta_info = basic_dict(centers, all_neighbors)
    test_crater_elevation(meta_info, xlen, ylen)

    return meta_info


def basic_dict(centers, neighbors):
    meta_info = {}
    c = 0
    m = len(centers)
    while c < m:
        meta_info[c] = {}
        meta_info[c]["neighbors"] = neighbors[c]
        meta_info[c]["center"] = centers[c]
        # since everything is square, if we have less than 4 neighboring cells
        # it's on the edge of the map
        if len(neighbors[c]) != 4:
            meta_info[c]["border"] = True
        else:
            meta_info[c]["border"] = False
        c += 1

    return meta_info


def recursive_drain(meta_info, river_tree_system, cell_list):
    for x in cell_list:
        if meta_info[x]["pool"] != None:
            meta_info[x]["pool"].riverdrain = True
        meta_info[x]["waterlevel"] = 0
        meta_info[x]["drained river"] = True

        if river_tree_system[x] == None:
            continue
        recursive_drain(meta_info, river_tree_system[x], river_tree_system[x])


def test_crater_elevation(meta_info, xlen, ylen):
    crater_center = np.array((7.5, 7.5, 0))
    radius = 5

    for c in meta_info:
        center = xlen/2, ylen/2
        p1 = np.array(meta_info[c]["center"])
        d = crater_center - p1

        my_mag = np.linalg.norm(d)
        elevation = 1-(abs(my_mag - 5)/radius)

        if my_mag < radius:
            elevation += 0.2
        if p1[0] > 7.5 and (6 < p1[1] < 9):
            elevation *= 0.6
        if elevation < 0.1:
            elevation = 0.1 + random.random() * 0.15

        meta_info[c]["elevation"] = elevation


def recalculate_slopes(meta_info):

    for cell in meta_info:
        slopes = {}
        biggest_diff = -float("inf")
        el1 = meta_info[cell]["elevation"]
        wl1 = meta_info[cell]["waterlevel"]

        t1 = el1 + wl1
        # my_el = meta_info[cell]["elevation"] + meta_info[cell]["waterlevel"]
        index = None

        steepest_slope = -float("inf")
        steepest_slope_id = None

        for neighborid in meta_info[cell]["neighbors"]:
            # actually put in the distance, not just the height difference.
            c1 = meta_info[cell]["center"]
            c2 = meta_info[neighborid]["center"]

            c1 = np.array(c1)
            c2 = np.array(c2)

            dvec = c2 - c1
            dist = np.linalg.norm(dvec)

            el2 = meta_info[neighborid]["elevation"]
            wl2 = meta_info[neighborid]["waterlevel"]

            t2 = el2 + wl2
            diff = t1 - t2

            if diff > steepest_slope and diff > 0:
                steepest_slope = diff
                steepest_slope_id = neighborid

            slope = diff / dist

            # this part is not needed for the calculations I'm doing
            # exactly here, but I think they're still useful
            # for a volume based approach.
            if diff > 0:
                # slope goes down hill.

                # static relevant difference
                H_sd = diff  # min(diff,wl1)

                # inertia resisting difference
                H_ir = 0

                # H_sr is pure speed relevant h
                if el1 > el2:
                    H_sr = wl1
                else:
                    H_sr = wl1 - H_ir

                slopes[neighborid] = (
                    diff, dvec, dist, slope, H_sd, H_sr, H_ir)
            else:
                # Don't do anything in particular,
                # this will be covered by this same function
                # but from the perspective of the other cell.

                # static relevant difference
                H_sd = 0
                # there is 0 static hydro based acceleration
                # from this cell / slope / neighbor

                # inertia resisting difference
                # this is the difference in terrain I'm flowing uphill
                # against
                H_ir = (el2-el1)

                # H_sr is pure speed relevant h
                if el1 > el2:
                    H_sr = wl1
                else:
                    H_sr = wl1 - H_ir

                slopes[neighborid] = (
                    diff, dvec, dist, slope, H_sd, H_sr, H_ir)

        meta_info[cell]["slopes"] = slopes
        if steepest_slope_id == None:
            meta_info[cell]["steepest slope"] = (None, None)
        else:
            meta_info[cell]["steepest slope"] = (
                steepest_slope_id, steepest_slope)


def rain_value_zeroing_step(meta_info, pools, rain_value=0.05):
    for cell in meta_info:
        old_water = float(meta_info[cell]["waterlevel"])
        meta_info[cell]["waterlevel"] *= 0.98  # evaporate some water
        meta_info[cell]["waterlevel"] += rain_value

        meta_info[cell]["waterdiff"] = 0
        meta_info[cell]["erosiondiff"] = 0
        meta_info[cell]["volume flow out"] = {}
        meta_info[cell]["volume flow in"] = {}
        meta_info[cell]["outflow point"] = False

        meta_info[cell]["drained river"] = False

        if cell in pools:
            meta_info[cell]["pool"] = pools[cell]
        else:
            meta_info[cell]["pool"] = None


def build_river_tree_map(meta_info):
    my_dict = {}

    temp = {}
    flow_direction = {}
    for x in meta_info:
        steepestid, steepest = meta_info[x]["steepest slope"]
        if steepestid != None:
            if steepestid not in flow_direction:
                flow_direction[steepestid] = [x]
            else:
                flow_direction[steepestid].append(x)

    visited = []
    full = {}
    for x in flow_direction:
        this = make_nested(flow_direction, x, visited)
        full[x] = this

    for x in visited:
        if x in full:
            full.pop(x)

    my_value_dict = {}
    depth_first_recursive_value_add(full, my_value_dict, 0)

    return full, my_value_dict, flow_direction


def make_nested(a, value, visited):
    full_ret_d = {}
    if value in a:
        for child in a[value]:
            ret_d = make_nested(a, child, visited)
            full_ret_d[child] = ret_d
            visited.append(child)
    else:
        full_ret_d = None

    return full_ret_d


def diff_apply(meta_info):
    for x in meta_info:
        meta_info[x]["waterlevel"] += meta_info[x]["waterdiff"]
        meta_info[x]["waterlevel"] = max(0, meta_info[x]["waterlevel"])
        meta_info[x]["elevation"] += meta_info[x]["erosiondiff"]
        meta_info[x]["elevation"] = max(0, meta_info[x]["elevation"])
        if meta_info[x]["border"]:
            meta_info[x]["waterlevel"] = 0


def calculate_diffs(meta_info, pools):
    for x in meta_info:
        if meta_info[x]["steepest slope"] == (None, None):
            if x not in pools:
                pools[x] = Pool(meta_info[x]["elevation"] +
                                meta_info[x]["waterlevel"])
                pools[x].member_cells.append(x)

            # nothing is flowing here.
            continue

        else:
            other_id, my_diff = meta_info[x]["steepest slope"]

        # abort and don't do anything if both are part of a pool.
        if x in pools and other_id in pools:
            # if it's the same one,
            if pools[x] == pools[other_id]:
                continue
            else:
                # if it's not the same one, combine them.
                this_pool = pools[x]
                for member in this_pool.member_cells:
                    pools[member] = pools[other_id]
                    meta_info[member]["pool"] = pools[other_id]
                continue

        # this is an outflow point
        if x in pools and other_id not in pools:
            # meta_info[x]["erosiondiff"] -= 0.005
            # meta_info[other_id]["erosiondiff"] += 0.005
            meta_info[x]["outflow point"] = True

        else:
            meta_info[x]["outflow point"] = False

        meta_info[x]["waterdiff"] -= my_diff * 0.2
        meta_info[other_id]["waterdiff"] += my_diff * 0.2


def depth_first_recursive_value_add(my_dict, value_dict, mydepth):

    value = 0
    for x in my_dict:

        if my_dict[x] == None:
            xvalue = 1
        else:
            xvalue = depth_first_recursive_value_add(
                my_dict[x], value_dict, mydepth+1)

        value_dict[x] = xvalue
        value += xvalue

    return value


def pool_integrity_check(meta_info, pools):
    # do I check everything?
    # yeah.

    keys = pools.keys()
    pool_list = []
    rml = []
    for x in keys:
        if pools[x] not in pool_list:
            if len(pools[x].member_cells) > 0:
                pool_list.append(pools[x])

    for pool in pool_list:
        # pick a random starting point, get all neighbors
        # if the neighborlist is shorter than my member_list, the difference is wrong.

        neighborlist = []

        cell = pool.member_cells[0]
        new = [cell]

        while True:
            new_new = []

            for x in new:
                for x2 in meta_info[x]["neighbors"]:
                    if x2 not in neighborlist and meta_info[x2]["waterlevel"] > 0:
                        new_new.append(x2)
                        neighborlist.append(x2)

            neighborlist += new
            new = new_new
            if new_new == []:
                break
        org_set = set(pool.member_cells)
        n_set = set(neighborlist)
        diff_list = list(org_set.difference(n_set))

        if diff_list != []:
            for x in diff_list:
                pools.pop(x)
                pool.member_cells.remove(x)

            this_pool = Pool(0)
            this_pool.member_cells = diff_list

            for x in this_pool.member_cells:
                pools[x] = this_pool
                meta_info[x]["pool"] = this_pool


def river_system_drain(meta_info, river_tree_system):
    for x in river_tree_system:
        # if it doesn't end at the end of the map,
        # just don't bother
        if meta_info[x]["border"] == False or meta_info[x]["pool"] != None:
            continue

        recursive_drain(meta_info, river_tree_system, [x])


def pool_loop(meta_info, pools):

    # so if I'm removing parts of my pool, I need to check when to split

    rml = []
    done = []
    pool_keys = list(pools.keys())
    for x in pool_keys:
        if x in done:
            continue
        my_pool = pools[x]
        smallrml = my_pool.equalize_main(meta_info)

        for x3 in meta_info:
            if meta_info[x3]["waterlevel"] <= 0:
                smallrml.append(x3)

        for x2 in smallrml:
            if x2 in pools:
                if x2 in pools[x2].member_cells:
                    pools[x2].member_cells.remove(x2)
                pools.pop(x2)
            done.append(x2)
        done.append(x)

    for x in pools:
        if x not in pools[x].member_cells:
            pools[x].member_cells.append(x)

    pool_integrity_check(meta_info, pools)


def clean_up_paths(my_paths, pools):
    new_paths = []
    for x in my_paths:
        new_path = []
        bool_list = []
        start_saving = False
        for x2 in x:
            bool_list.append(x2 in pools)

        assert len(bool_list) == len(x)
        c = 0
        while c < len(bool_list)-1:
            # print(x[c],bool_list[c])
            if start_saving:
                new_path.append(x[c])
            if bool_list[c] == True and bool_list[c+1] == False:
                new_path += [x[c], x[c+1]]
                start_saving = True
                c += 1

            elif bool_list[c] == False and bool_list[c+1] == True:
                start_saving = False
                new_path.append(x[c+1])
                c += 1
            c += 1
        if new_path != [] and len(new_path) > 1:
            new_paths.append(new_path)
    # new_new_paths = []
    # for path in new_paths:
     #   new_new_paths.append(tuple(path))
    # new_paths= list(set(new_new_paths))
    return new_paths


def find_lake_draining_paths(meta_info, pools, local_river_tree_level, current_level):
    """this finds rivers that go lake to lake."""
    my_completed_lists = []

    for x in local_river_tree_level:

        # if local_river_tree_level[x] == None:
        # this is an end, and potentially something that I want.
        if meta_info[x]["pool"] != None or x in pools:
            this = list(current_level)
            this.append(x)

            if len(this) > 1:
                my_completed_lists.append(this)

        if local_river_tree_level[x] != None:
            new_branch_list = list(current_level)
            new_branch_list.append(x)

            r = find_lake_draining_paths(
                meta_info, pools, local_river_tree_level[x], new_branch_list)

            my_completed_lists += r

    return my_completed_lists


def find_lake_to_edge_draining_paths(meta_info, pools, local_river_tree_level, current_level):
    """this finds rivers that go lake to lake."""
    my_completed_lists = []

    for x in local_river_tree_level:

        this = list(current_level)
        this.append(x)
        is_pool = (meta_info[x]["pool"] != None or x in pools)
        is_end = local_river_tree_level[x] == None

        if not is_end:
            new_branch_list = list(current_level)
            new_branch_list.append(x)

            r = find_lake_draining_paths(
                meta_info, pools, local_river_tree_level[x], new_branch_list)

            my_completed_lists += r

        if len(this) > 1 and is_pool and is_end:
            my_completed_lists.append(this)

    return my_completed_lists


def river_gradient_erosion(meta_info, my_paths):
    for path in my_paths:
        els = []
        wls = []
        pathlen = len(path)
        if pathlen > 2:
            wl_sum = 0
            for index in path:
                el = meta_info[index]["elevation"]
                wl = meta_info[index]["waterlevel"]
                els.append(el)
                wls.append(wl)
            wl_sum = sum(wls)
            mm = els[-1]
            mn = els[0]

            overall_gradient_step = (mm-mn) / (pathlen-1)
            c = 0
            while c < pathlen:
                index = path[c]
                el = meta_info[index]["elevation"]
                meta_info[index]["waterlevel"] = wl_sum / pathlen

                goal = (mn+c*overall_gradient_step)
                diff = goal - el

                # if the elevation is smaller than the goal, add soemthing
                # else the difference is negative and the whole thing will be
                # lowered by itself.
                
                meta_info[index]["elevation"] += diff * 0.6

                c += 1


def set_river_flag(meta_info,my_paths):
    for x in meta_info:
        meta_info[x]["part of a river"]=False
        for path in my_paths:
            if x in path:
                meta_info[x]["part of a river"]=True

def clear_pool_members_via_paths(meta_info, path_list):
    """this is for multiple lists of paths"""
    # my_paths
    for path_set in path_list:
        for path in path_set:
            for x in path:
                this_pool=meta_info[x]["pool"]
                if this_pool!=None:
                    if x in this_pool.member_cells:
                        this_pool.member_cells.remove(x)
                meta_info[x]["pool"]=None
