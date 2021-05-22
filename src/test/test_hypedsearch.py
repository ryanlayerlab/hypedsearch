import hypedsearch

def test_add_pass():
    assert 1+1 == 2, "Should be 2"

def test_add_fail():
    assert 1+1 == 1, "Should be 2"

if __name__ == "__main__":
    test_add_pass()
    test_add_fail()
    print("All tests passed")
