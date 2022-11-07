from pyoptools.raytrace.library import library
from pyoptools.raytrace._comp_lib.optic_factory import optic_factory

def main():

    failed = 0
    for part, descriptor in library.items():
        #print('Part : ', part)
        try:
            optic = optic_factory(**descriptor)
        except Exception:
            print('Failed : ', part)
            failed += 1

    print('failed : ', failed)

if __name__ == '__main__':
    main()
