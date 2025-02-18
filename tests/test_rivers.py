import unittest
from simple_rivers import simple_rivers

class TestMyRivers(unittest.TestCase):
    
    def test_basic(self):
        simple_rivers.full_simple_river_main()
    
    def test_find_gradients(self):
        
        elevation = {1:5,2:10}
        waterlevels = {1:0,2:0}
        neighbors = {1:[2],2:[1]}
        centers = {1:(0,0,0),2:(1,0,0)}
        
        
        slopes, steepest_slopes = simple_rivers.find_gradients(elevation,waterlevels,neighbors,centers)
        
        assert slopes == {1: {2: (-5, (1, 0, 0), 1.0, -5.0, 0, 5, 5)}, 2: {1: (5, (-1, 0, 0), 1.0, 5.0, 5, 0, 0)}}
        assert steepest_slopes == {1: (2, -5), 2: (1, 5)}
    
    def test_build_river_tree_map(self):
        # this tells me nothing... meh.
        
        meta_info = {}
        r = simple_rivers.build_river_tree_map(meta_info)
    
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
