import unittest
from unittest import TestCase
from circuit_benchmarks.qft import *

class TestQFT(TestCase):
    @unittest.skip("not a test")
    def test_print_qft(self) -> None:
        c = qft(3)
        print(c)

