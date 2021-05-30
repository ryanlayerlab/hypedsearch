import sys
import argparse

def string_to_bool(s: str) -> bool:
    s = str(s)
    if s.lower() == 'false' or 'f' in s.lower():
        return False
    return True

def main():
    print("Hello World!")

if __name__ == "__main__":
    main()
