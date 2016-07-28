#!/bin/bash

g++ -std=c++11 -o flow_tree flow_tree.C -lstdc++ `root-config --libs` -I$ROOTSYS/include

g++ -std=c++11 -o analyze analyze.C -lstdc++ `root-config --libs` -I$ROOTSYS/include


