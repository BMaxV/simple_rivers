
import uuid

from my_save import sxml_main
from vector import vector

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

class dummypoly:
    def __init__(self):
        self.center = vector.Vector(0,0,0)
        self.points = []

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
            neighbors = {}
            for p in directions:
                if p in centers:
                    my_index = centers.index(p)
                    neighbors[my_index] = []
            all_neighbors[total_counter] = neighbors
            total_counter += 1
            y += 1
        x += 1
        
    meta_info = basic_dict(centers, all_neighbors)
    test_crater_elevation(meta_info, xlen, ylen)

    return meta_info

def find_starting_points(river_trees):
    targets = []
    keys = list(river_trees.keys())
    for x in river_trees:
        targets += list(river_trees[x])
    
    targets = set(targets)
    keys = set(keys)
    
    print(len(targets))
    print(len(keys))
    
    # so these are the keys that are not being targeted.
    starting_points = keys.difference(targets)
    starting_points = list(starting_points)
    return starting_points

def basic_dict(centers, neighbors):
    meta_info = {}
    c = 0
    m = len(centers)
    while c < m:
        meta_info[c] = {}
        
        meta_info[c]["polygon"] = dummypoly()
        
        this_center = vector.Vector(*centers[c])
        points = []
        diffs = [vector.Vector(1,1,0)*0.5,
                vector.Vector(1,-1,0)*0.5,
                vector.Vector(-1,-1,0)*0.5,
                vector.Vector(-1,1,0)*0.5]
                
        for x in diffs:
            points.append(tuple(this_center+x))
            
        meta_info[c]["polygon"].points = points
        meta_info[c]["waterlevel"] = 0
        meta_info[c]["neighbors"] = neighbors[c]
        meta_info[c]["center"] = centers[c]
        # since everything is square, if we have less than 4 neighboring cells
        # it's on the edge of the map
        
        
        
        if len(neighbors[c]) != 4:
            meta_info[c]["border"] = True
        else:
            meta_info[c]["border"] = False
        
        c += 1
    
    for cell_index in meta_info:
        points1 = meta_info[cell_index]["polygon"].points
        for ni in meta_info[cell_index]["neighbors"]:
            points2 = meta_info[ni]["polygon"].points
            c2 = -1
            list_of_tuples = []
            while c2 < len(points1)-1:
                pi1 = points1[c2]
                pi2 = points1[c2+1]
                if pi1 in points2 and pi2 in points2:
                    # this avoids -1 stuff
                    index1 = points1.index(pi1)
                    index2 = points1.index(pi2)
                    list_of_tuples.append((index1,index2))
                c2 += 1
            meta_info[cell_index]["neighbors"][ni] = list_of_tuples
    
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
    crater_center = vector.Vector(7.5, 7.5, 0)
    radius = 5

    for c in meta_info:
        center = xlen/2, ylen/2
        p1 = vector.Vector(*meta_info[c]["center"])
        d = crater_center - p1

        my_mag = d.magnitude()
        elevation = 1-(abs(my_mag - 5)/radius)

        if my_mag < radius:
            elevation += 0.2
        if p1[0] > 7.5 and (6 < p1[1] < 9):
            elevation *= 0.6
        if elevation < 0.1:
            elevation = 0.1

        meta_info[c]["elevation"] = elevation


def unpack_meta(meta_info):
    elevation = {}
    waterlevel = {}
    neighbors = {}
    centers = {}
    for x in meta_info:
        elevation[x]=meta_info[x]["elevation"]
        waterlevel[x]=meta_info[x]["waterlevel"]
        neighbors[x]=meta_info[x]["neighbors"]
        centers[x]=meta_info[x]["center"]
    return elevation, waterlevel, neighbors, centers

def merge_slopes(meta_info,slopes,steepest_slopes):
    for key in slopes:
        meta_info[key]["slopes"]=slopes[key]
        meta_info[key]["steepest slope"]=steepest_slopes[key]

def recalculate_slopes(meta_info):
    """
    this actually doesn't take too many inputs,
    they are just spread around in an unfortunate way.
    """
    
    elevation, waterlevel, neighbors, centers = unpack_meta(meta_info)
    slopes, steepest_slopes = find_gradients(elevation, waterlevel, neighbors, centers)
    merge_slopes(meta_info,slopes,steepest_slopes)
    
def find_gradients(elevation, waterlevel, neighbors, centers):
    slopes = {}
    steepest_slopes = {}
    for cell in elevation:
        small_slopes = {}
        
        el1 = elevation[cell]
        wl1 = waterlevel[cell]

        t1 = el1 + wl1
        
        index = None

        steepest_slope = -float("inf")
        steepest_slope_id = None
        slope_list = []
        # I need this to exist.
        # I think. 90% sure.
        assert len(neighbors[cell]) > 0
        
        for neighborid in neighbors[cell]:
            # actually put in the distance, not just the height difference.
            
            c1 = centers[cell]
            c2 = centers[neighborid]

            c1 = vector.Vector(*c1)
            c2 = vector.Vector(*c2)

            dvec = c2 - c1
            dist = dvec.magnitude()

            el2 = elevation[neighborid]
            wl2 = waterlevel[neighborid]

            t2 = el2 + wl2
            diff = t1 - t2

            slope_list.append((neighborid,diff))

            slope = diff / dist

            # this part is not needed for the calculations I'm doing
            # exactly here, but I think they're still useful
            # for a volume based approach.
            
            same_tup = [diff, dvec, dist, slope]
            h_stuff = get_H_stuff(diff, el1, el2, wl1, wl2)
            
            slop_tup = tuple(same_tup + h_stuff)
            small_slopes[neighborid] = slop_tup
        
        steepest_slope_id, steepest_slope = max(slope_list,key=lambda x:x[1])
        
        slopes[cell] = small_slopes
        if steepest_slope_id == None:
            steepest_slopes[cell] = (None, None)
        else:
            slop_tup2 = (steepest_slope_id, steepest_slope)
            steepest_slopes[cell] = slop_tup2
    
    return slopes, steepest_slopes
    

def get_H_stuff(diff, el1, el2, wl1, wl2):
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
    h_stuff = [H_sd, H_ir, H_ir]
    return h_stuff

def full_simple_river_main():
    
    meta_info = generate_basic_data()
    pools = {}
    all_are_border_tiles = False
    time_c = 0
    my_max = 60
    while all_are_border_tiles == False and time_c <  my_max:
        all_are_border_tiles, my_trees = main_step(meta_info,pools,all_are_border_tiles)
        
        # if you want to do something each step, e.g.
        # producing some kind of step by step output, do that here.
        
        time_c += 1
        
    #output_rivers(meta_info,my_trees,time_c)
    make_save(meta_info,my_trees)
    return meta_info,my_trees



def simple_river_main(meta_info):
    
    pools = {}
    all_are_border_tiles = False
    time_c = 0
    my_max = 60
    while all_are_border_tiles == False and time_c <  my_max:
        #print("yep")
        all_are_border_tiles, my_trees = main_step(meta_info,pools,all_are_border_tiles)
        
        # if you want to do something each step, e.g.
        # producing some kind of step by step output, do that here.
        
        time_c += 1
        
    print("stopped after",time_c,"out of",my_max)
        
    #output_rivers(meta_info,my_trees,time_c)
    #make_save(meta_info,my_trees)
    return meta_info,my_trees

def make_save(meta_info,my_trees):
    save_dict = {}
    for cell in meta_info:
        save_dict[cell] = {}
        l = []
        l2 = meta_info[cell]["polygon"].points
        
        for x in l2:
            v = []
            for x2 in x:
                v.append(float(x2))
            l.append(tuple(v)) 
        
        save_dict[cell]["points"] = l
        save_dict[cell]["neighbors"]  = meta_info[cell]["neighbors"]
        save_dict[cell]["elevation"] = float(meta_info[cell]["elevation"])
        save_dict[cell]["river value"]  = meta_info[cell]["river value"]
                
    save_dict["rivertrees"] = my_trees
    fn = "saved_rivers.xml"
    save_dict = {"data":save_dict}
    sxml_main.write(fn,save_dict)

def main_step(meta_info,pools,all_are_border_tiles):
    """
    this main step is meant to be run inside of a loop
    if "all_are_border_tiles" is true, "my trees" will contain
    a full gradient tree for every single tile in map,
    describing how all water will drain.
    
    if there are depressions that would fill  with water, I will
    fill them with water, determine the outflow point and
    erode the outflow path until the lowest point of that lake can drain
    to the edge of the map
    
    if you want to keep a few lakes,
    the best approach is probably to let this function run,
    note the number of steps it took, and then run it again and stop early
    
    If you have specified a seed and the whole process is deterministic,
    as it should be, there should be some lakes in the data.
    """
    # how will things flow?
    
    
    
    #recalculate_slopes(meta_info)
    
    elevation, waterlevel, neighbors, centers = unpack_meta(meta_info)
    slopes, steepest_slopes = find_gradients(elevation, waterlevel, neighbors, centers)
    merge_slopes(meta_info,slopes,steepest_slopes)
    
    my_trees, depth_data, base_dict = build_river_tree_map(steepest_slopes)
    
    # do all tiles drain to the edge of the map somehow?
    all_are_border_tiles = True
    for x in my_trees:
        if meta_info[x]["border"] == False:
            all_are_border_tiles = False
            break
    
    # zero all diffs for this, apply some rain
    rain_value_zeroing_step(meta_info,pools,0.01)
    
    for x in meta_info:
        if x in depth_data:
            meta_info[x]["river value"] = depth_data[x]
        else:
            meta_info[x]["river value"] = 1 # weird?
            
    # equalize the water in lakes / pools,#
    # calculate flow diffs
    # apply them
    # equalize again
    pool_loop(meta_info,pools)
    
    calculate_diffs(meta_info,pools)
    diff_apply(meta_info)
    
    pool_loop(meta_info,pools)
    
    # maybe I have to reroute my drainage through my pools.
    
    # drain all water from all regular rivers and tiles
    river_system_drain(meta_info,my_trees)
    
    # when I'm doing this and the last cell of that path is a 
    # lake cell, I'm draining the lake cell, but I'm not cutting down
    # the path.
    
    # also, it might be an idea to include lakes into the 
    # river draining system.
    
    # go through my river / lake system and also drain
    # all lakes that have some connection to the edge of the map
    
    current_level = []
    my_paths = find_lake_draining_paths(meta_info,pools,my_trees,current_level)
    my_paths = clean_up_paths(my_paths,pools)
    
    current_level = []
    lake_edge_drain_paths = find_lake_to_edge_draining_paths(meta_info,pools,my_trees,current_level)
    
    # path erosion
    # this will erode the river paths
    # 
    river_gradient_erosion(meta_info,my_paths)
    river_gradient_erosion(meta_info,lake_edge_drain_paths)
    path_list = [my_paths,lake_edge_drain_paths]
    
    clear_pool_members_via_paths(meta_info,path_list)

    # I can also just not add rain water to cells that are being drained.
    set_river_flag(meta_info,my_paths)
    return all_are_border_tiles, my_trees

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

def find_flow_direction(steepest_slopes):
    
    flow_direction = {}
    for cell_id in steepest_slopes:
        steepestid, steepest = steepest_slopes[cell_id]
        if steepestid != None:
            if steepestid not in flow_direction:
                flow_direction[steepestid] = [cell_id]
            else:
                flow_direction[steepestid].append(cell_id)
    return flow_direction

def make_nested_new(layer,inverted_flow):
    sub_tree = {}
    for x in layer:
        if x in inverted_flow:
            new_layer=inverted_flow[x]
            r=make_nested_new(new_layer,inverted_flow)
            sub_tree[x]=r
        else:
            sub_tree[x]=None
    return sub_tree

def build_river_tree_map(steepest_slopes):
    
    # so this is a dict, pointing to list of ids, where
    # it's flowing. but I think... these lists are just 1 long?
    #flow_direction = find_flow_direction(steepest_slopes)
    
    inverted_flow = {}
    for key1 in steepest_slopes:
        other_id, d_h = steepest_slopes[key1]
        inverted_flow[other_id] = key1
    
    keys1 = set(list(steepest_slopes.keys()))
    keys2 = set(list(inverted_flow.keys()))
    roots = list(keys1.difference(keys2))
    
    layer = list(roots)
    trees = make_nested_new(layer,inverted_flow)
        # that's dumb, I don't want to remove, I want to never add this.
    
    my_value_dict = {}
    depth_first_recursive_value_add(trees, my_value_dict, 0)

    return trees, my_value_dict, inverted_flow


def make_nested(flow_direction, value, visited):
    full_ret_d = {}
    if value in flow_direction:
        for child in flow_direction[value]:
            ret_d = make_nested(flow_direction, child, visited)
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
    """what does this do?"""
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
        
        c = 0
        m = 1000
        while True and c < 1000:
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
            c+=1
        
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
    """
    split paths based on whether they go through pools or not.
    """
    new_paths = []
    for path in my_paths:
        
        new_path = []
        bool_list = []
        save_section = False
        
        for cell in path:
            bool_list.append(cell in pools)

        c = 0
        while c < len(bool_list)-1:
            
            this_in_pool = (bool_list[c] == True)
            next_not_in_pool = (bool_list[c+1] == False)
            
            if save_section:
                new_path.append(path[c])
            
            # turn it on
            if this_in_pool and next_not_in_pool:
                new_path += [path[c], path[c+1]]
                save_section = True
                c += 1
            
            # turn it off
            elif not this_in_pool and not next_not_in_pool:
                save_section = False
                new_path.append(path[c+1])
                c += 1
                
            c += 1
            
        if new_path != [] and len(new_path) > 1:
            new_paths.append(new_path)
    
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

def unpack_river_values_neighbors(d):
    neighbor_dict = {}
    river_value_dict = {}
    for cell in d:
        sub_dict={}
        for x in d[cell]["neighbors"]:
            neighbor_value = d[cell]["neighbors"][x]
            sub_dict[x]= neighbor_value[0]
        neighbor_dict[cell] = sub_dict
        
        river_value_dict[cell]=d[cell]["river value"]
    return neighbor_dict, river_value_dict

def build_my_edge_list(neighbor_dict, river_value_dict, index, my_river_tuple_list,threshold_value=1):
    """
    I have code that's creating custom objects based
    on the data relative to the neighbor, e.g. 
    
    "
    do I want a door between tile A and tile B? 
    If (value) > something: yes. 
    Then record the edge between A and B.
    "
    
    flat is the bool for "no, don't do anything"
    
    really this needs a river value dict and a neighbor dict.
    
    
    hmmmmmmmmmmm.
    the mydictionary [cell][propertyname1] -> value1
    the mydictionary [cell][propertyname2] -> value2
    is really convenient for the rivers,
    but less so for the planet, where I'm layering 
    different systems.
    
    
    
    
    """
    
    flat = True
    my_edge_list = []
    if river_value_dict[index] > threshold_value:
        flat = False
        my_edge_list = []
        for n_index in neighbor_dict[index]:
            if (index,n_index) not in my_river_tuple_list:
                continue
            
            edge_tuple = neighbor_dict[index][n_index]
            
            if river_value_dict[n_index] > threshold_value:
                my_edge_list += [list(edge_tuple)]
                
    return my_edge_list, flat

def recursive_river_tree_unpack(my_dict):
    
    my_tuple_list = []
    for x in my_dict:
        if my_dict[x]!=None:
            lower_list = recursive_river_tree_unpack(my_dict[x])
            other_keys = list(my_dict[x].keys())
            for x2 in other_keys:
                my_tuple_list.append((x,x2))
                my_tuple_list.append((x2,x))
            my_tuple_list += lower_list
            
    return my_tuple_list


if __name__=="__main__":
    full_simple_river_main()
