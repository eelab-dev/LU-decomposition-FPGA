# Clean generated programmes without file extension

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
path = os.getcwd()
for files in [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]:
    temp = os.path.splitext(files)
    if not temp[1]:
        os.remove(temp[0])
        print(f"Delete '{temp[0]}'")
