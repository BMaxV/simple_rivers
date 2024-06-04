import os
import setuptools

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setuptools.setup(
    name = "simple_rivers",
    version = "0.10",
    author = "Brumo Maximilian Voss",
    author_email = "bruno.m.voss@gmail.com",
    description = ("carve some simple rivers into your proc gen landscape"),
    packages=["simple_rivers"],
    license = "MIT",
    #keywords = "sympy 3d",
   
)
