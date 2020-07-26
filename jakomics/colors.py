# https://stackoverflow.com/questions/287871/how-to-print-colored-text-in-terminal-in-python


def test():
    print("colors module loaded correctly")


class bcolors:
    PURPLE = '\033[95m'
    CYAN = "\033[36m"
    BLUE = '\033[94m'
    YELLOW = '\033[93m'
    GREEN = '\033[92m'
    RED = '\033[91m'
    END = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
