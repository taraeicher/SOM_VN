import pysam
import sys

def main():
    file_name = sys.argv[1]
    print(file_name)
    pysam.index(file_name)
    
if __name__ == "__main__":
    main()