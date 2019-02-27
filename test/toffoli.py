import unittest
from unittest import TestCase
from circuit_benchmarks.toffoli import *

class TestToffoli(TestCase):
    # @unittest.skip("not a test")
    def test_print_toffoli(self) -> None:
        c = toffoli(2)
        print(c)

