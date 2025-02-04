from src.erik_utils import run_in_parallel


class TestRunInParallel:
    @staticmethod
    def test_add_five():
        fcn = lambda x: x + 5
        input_array = [1, 2, 3, 4, 5]
        actual = list(run_in_parallel(fcn_of_one_variable=fcn, input_array=input_array))
        assert actual == [fcn(x) for x in input_array]
