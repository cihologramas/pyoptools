import os


def findfiles2clean(dir,files=[], ):
    for file in os.listdir(dir):
        if file != "build":
            path = os.path.join(dir, file)
            if os.path.isfile(path) and path.endswith(".pyx"):
                files.append(path[:-4]+".c")
                files.append(path[:-4]+".so")
            elif os.path.isdir(path):
                findfiles2clean(path, files)
    return files


f2c = findfiles2clean(".")
for i in f2c:
    print("Deleting: ", i)
    try:
        os.remove(i)
    except:
        pass
