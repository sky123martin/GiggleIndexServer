import glob
import os
from config import config

config = config

def delete_directory(dir):  # Recursively deletes directory
    if os.path.isdir(dir):
        print("Deleting", dir)
        item_list = glob.glob(os.path.join(dir, "*"))
        for item in item_list:
            if os.path.isfile(item):
                os.remove(item)
            elif os.path.isdir(item):
                delete_directory(item) 
        os.rmdir(dir)

if __name__ == "__main__":
    delete_directory("data/")
    delete_directory("indices/")
    if os.path.isfile(config.DB):
        os.remove(config.DB)