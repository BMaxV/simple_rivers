import unittest
from simple_rivers import simple_rivers

class TestMyRivers(unittest.TestCase):
    
    def test_basic(self):
        simple_rivers.full_simple_river_main()
    
    
    def test_get_roots_branches(self):
        
        up = {1: [2, 3, 4], 4: [5]}
        down = {2: 1, 3: 1, 4: 1, 5: 4}
        roots, branches = simple_rivers.get_roots_branches(up,down)
        
        assert roots == [1]
        assert branches == [2,3,5]
        
    def test_find_gradients(self):
        
        elevation = {1:5,2:10,3:7,4:7,5:8}
        waterlevels = {1:0,2:0,3:0,4:0,5:0}
        neighbors = {1:[2,3,4],2:[1,3,4],4:[1,2,5],3:[1,2],5:[4]}
        centers = {1:(0,0,0),2:(1,0,0),3:(0.5,0.5,0),4:(0.5,-0.5,0),5:(-0.6,-0.6,0)}
        
        
        slopes, steepest_slopes, up, down = simple_rivers.find_gradients(elevation,waterlevels,neighbors,centers)
        
        assert slopes == {1: {2: (-5, (1, 0, 0), 1.0, -5.0, 0, 5, 5), 3: (-2, (0.5, 0.5, 0), 0.7071067811865476, -2.82842712474619, 0, 2, 2), 4: (-2, (0.5, -0.5, 0), 0.7071067811865476, -2.82842712474619, 0, 2, 2)}, 2: {1: (5, (-1, 0, 0), 1.0, 5.0, 5, 0, 0), 3: (3, (-0.5, 0.5, 0), 0.7071067811865476, 4.242640687119285, 3, 0, 0), 4: (3, (-0.5, -0.5, 0), 0.7071067811865476, 4.242640687119285, 3, 0, 0)}, 3: {1: (2, (-0.5, -0.5, 0), 0.7071067811865476, 2.82842712474619, 2, 0, 0), 2: (-3, (0.5, -0.5, 0), 0.7071067811865476, -4.242640687119285, 0, 3, 3)}, 4: {1: (2, (-0.5, 0.5, 0), 0.7071067811865476, 2.82842712474619, 2, 0, 0), 2: (-3, (0.5, 0.5, 0), 0.7071067811865476, -4.242640687119285, 0, 3, 3), 5: (-1, (-1.1, -0.1, 0), 1.104536101718726, -0.9053574604251853, 0, 1, 1)}, 5: {4: (1, (1.1, 0.1, 0), 1.104536101718726, 0.9053574604251853, 1, 0, 0)}}
        assert steepest_slopes == {1: (None, None), 2: (1, 5), 3: (1, 2), 4: (1, 2), 5: (4, 1)}

        roots, branches = simple_rivers.get_roots_branches(up,down)
    
    def test_get_roots_branches(self):
        
        elevation = {1:5,2:10,3:7,4:7,5:8}
        waterlevels = {1:0,2:0,3:0,4:0,5:0}
        neighbors = {1:[2,3,4],2:[1,3,4],4:[1,2,5],3:[1,2],5:[4]}
        centers = {1:(0,0,0),2:(1,0,0),3:(0.5,0.5,0),4:(0.5,-0.5,0),5:(-0.6,-0.6,0)}
        
        slopes, steepest_slopes, up, down = simple_rivers.find_gradients(elevation,waterlevels,neighbors,centers)
        
        roots, branches = simple_rivers.get_roots_branches(up,down)
        
        assert roots == [1]
        assert branches == [2, 3, 5]


        
    def test_build_river_tree_map(self):
        
        elevation = {1:5,2:10,3:7,4:7,5:8}
        waterlevels = {1:0,2:0,3:0,4:0,5:0}
        neighbors = {1:[2,3,4],2:[1,3,4],4:[1,2,5],3:[1,2],5:[4]}
        centers = {1:(0,0,0),2:(1,0,0),3:(0.5,0.5,0),4:(0.5,-0.5,0),5:(-0.6,-0.6,0)}
        
        slopes, steepest_slopes, up, down = simple_rivers.find_gradients(elevation,waterlevels,neighbors,centers)
        
        roots, branches = simple_rivers.get_roots_branches(up,down)
        trees, my_value_dict = simple_rivers.build_river_tree_map(roots,up)
        
        assert my_value_dict =={2: 1, 3: 1, 5: 1, 4: 1, 1: 3}
        assert trees == {1: {2: None, 3: None, 4: {5: None}}}
        
    def test_make_nested(self):
        
        layer = [1]
        inverted = {1:[2],
                    2:[3],}
   
        trees = simple_rivers.make_nested_new(layer, inverted)
        assert trees == {1:{2:{3:None}}}
    
def test_single():
    single = TestMyRivers()
    
    #single.test_partial_generation_save_then_load()
    a=1

if __name__=="__main__":
    unittest.main()
