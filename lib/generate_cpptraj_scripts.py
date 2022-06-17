#!/usr/bin/env python

import numpy as np
import sys
import yaml

if __name__ == '__main__':
    if len(sys.argv) !=2:
        print('Usage: generate_cpptraj_scripts.py input.yaml')
        exit(1)
    yaml_input = sys.argv[1]
