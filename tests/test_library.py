from pyoptools.raytrace.library import library

def main():

    thor = library.Thorlabs

    #for k, v in library.items():
    #    if k == '08083' or k == 'LB1844-C':
    #        print(k, v['description'])

    #d = library['08083']

    lens = thor.get('LB1862')
    print(type(lens))
    print(lens)
    print('\n')

    lens2 = library.get('LB1862')
    print(type(lens2))
    print(lens2)

if __name__ == '__main__':
    main()
