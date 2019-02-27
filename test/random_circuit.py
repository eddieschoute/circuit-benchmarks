import unittest
from unittest import TestCase
from circuit_benchmarks.random_circuit import *

class TestRandomCircuit(TestCase):
    @unittest.skip("not a test")
    def test_print_random_circuit(self) -> None:
        c = random_circuit(10)
        print(c)

