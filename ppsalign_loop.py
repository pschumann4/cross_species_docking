####################################################################################################
# This function will take a file directory of .poc files as an input and run PPSalign on each.     #
# The user will need to set their current working directory to the directory containing the        #
# .poc files and have the PPSalign executable in their PATH prior to running this function.        #
####################################################################################################

import os

def ppsalign_loop():
    # Get the path to the directory containing the .poc files
    poc_dir = os.getcwd()
    poc_files = os.listdir(poc_dir)

    # Specify name of template .poc file
    template = input('Enter the name of the template .poc file: ')
    if not template.endswith('.poc'):
        template += '.poc'
    
    for poc in poc_files:
        if poc.endswith('.poc') and poc != template:
            # Create the command to run PPSalign
            cmd = 'PPSalign ' + poc + ' ' + template + ' > ' + poc.split('.poc')[0] + '_PPS.txt'
            os.system(cmd)
        else:
            continue
    
    # Move all of the _PPS.txt files to a folder called 'PPS_files'
    if not os.path.isdir('PPS_files'):
        os.mkdir('PPS_files')
    for file in os.listdir(poc_dir):
        if file.endswith('_PPS.txt'):
            os.rename(file, 'PPS_files/' + file)
    print('Done!')

if __name__ == '__main__':
    ppsalign_loop()