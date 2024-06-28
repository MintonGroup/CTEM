"""
 Copyright 2024 - The Minton Group at Purdue University
 This file is part of CTEM
 CTEM is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 CTEM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
 of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 You should have received a copy of the GNU General Public License along with CTEM. 
 If not, see: https://www.gnu.org/licenses. 
"""

import ctem
import unittest
import os
import numpy as np
import tempfile

class TestCTEMIO(unittest.TestCase):
    def setUp(self):
        # Initialize a target and surface for testing
        self.tmpdir=tempfile.TemporaryDirectory()
        self.simdir = self.tmpdir.name
        
    def tearDown(self):
        # Clean up temporary directory
        self.tmpdir.cleanup() 

    def test_single_crater(self):
        """
        Tests that CTEM is able to successfully generate a test crater without any exceptions being raised
        """
        
        sim = ctem.Simulation()
        try:
            sim.run()
        except Exception as e:
            self.fail(f"CTEM raised an exception: {e}")

        return
    
    
if __name__ == '__main__':
    unittest.main()