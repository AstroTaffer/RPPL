from os import mkdir

mkdir("backup/")
with open("backup/testfile.txt", "w") as testfile:
    testfile.write("test_string")


"""

"""
