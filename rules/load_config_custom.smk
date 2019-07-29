import os

configfile: "config/customize_data.yaml"


class PathNotGiven(RuntimeError):
    def __init__(self, message):
        self.message = message


try:
    refs = list(map(str.strip, config["refs"].split(",")))
    project_dir = config["outpath"].rstrip("/")
except AttributeError as e:
    #print("Please specify the reference genome files and output directory.", e)
    raise PathNotGiven(
        "The reference genome files or output directory are not specified.")

threads = config["threads"]
results_dir = project_dir + "/results"
