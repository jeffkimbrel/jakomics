import sys
import os
import argparse
import time
from jakomics import utilities, colors

print(f'{colors.bcolors.YELLOW}THIS CODE IS UNDER DEVELOPMENT!{colors.bcolors.END}')

# utilities.test()
# colors.test()

# OPTIONS #####################################################################

parser = argparse.ArgumentParser(description="", formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--in_dir',
                    help="Directory with faa files",
                    required=False,
                    default="")

parser.add_argument('-f', '--files',
                    help="Paths to individual faa files",
                    nargs='*',
                    required=False,
                    default=[])

parser.add_argument('--out_dir',
                    help="Directory to write results to",
                    required=True)

parser.add_argument('--threads',
                    help="Threads",
                    default=8,
                    type=int,
                    required=False)

args = parser.parse_args()

args.out_dir = os.path.abspath(args.out_dir) + '/'

if not os.path.exists(args.out_dir):
    print("\nCreating directory " + args.out_dir)
    os.makedirs(args.out_dir)

counter = 0

# FUNCTIONS ###################################################################


def main(input_path, output_path):
    global counter
    file_name = os.path.basename(input_path)
    output_file = os.path.join(output_path, file_name + '.kofam93.txt')
    output_temp = os.path.join(output_path, file_name + '.tmp')

    if os.path.exists(output_file):
        print(f'{colors.bcolors.RED}Skipping {output_file} because it already exists{colors.bcolors.END}')
    else:
        command = '/Users/kimbrel1/Science/kofam_93/kofamscan-1.1.0/exec_annotation -f mapper --tmp-dir ' + \
            output_temp + ' -o ' + output_file + ' ' + input_path + ' --cpu 8'

        print(
            f'Finished {counter} of {len(file_list)}... currently running {colors.bcolors.GREEN}{file_name}{colors.bcolors.END}\t\t\t',
            end="\r",
            file=sys.stderr)

        os.system(command)
        os.system('rm -fR ' + output_temp)
        time.sleep(2)  # sleep so control c will cancel easier

    counter += 1

## MAIN LOOP ###################################################################


if __name__ == "__main__":

    file_list = utilities.get_file_list(args.files, [''])
    file_list = utilities.get_directory_file_list(args.in_dir, [''], file_list)

    if len(file_list) == 0:
        sys.exit(
            f"{colors.bcolors.RED}Error: No valid .faa files were found!{colors.bcolors.END}")

    for file in file_list:
        main(file, args.out_dir)

    print()
