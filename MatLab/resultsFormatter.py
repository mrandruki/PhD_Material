from datetime import datetime
from imageio import imread
from numpy import zeros, array, uint8
from scipy.io import loadmat, savemat
from time import sleep

imagePath = input('Please enter the name of the image file to process: ')
imageName = imagePath[:-4]
imageExtension = imagePath[-4:]

## Andy's Laptop
matlabPath = r'C:\Users\Tesla\Dropbox\PhD\MatLab\\'

## Andy's PC
##matlabPath = r'C:\Users\Andy\Dropbox\PhD\MatLab\\'

## Drew's
##matlabPath = r'C:\Users\dange\Dropbox\PhD\MatLab\\'

data = {}
data['__header__'] = b'MATLAB 5.0 MAT-file, Platform: PCWIN64, Created on: ' + \
                     str.encode(datetime.now().strftime("%a %b %d %H:%M:%S %Y"))
data['__version__'] = '1.0'
data['__globals__'] = ['BOUNDARY', 'FLOW', 'SOLID']
data['BOUNDARY'] = array([[0]], dtype=uint8)
data['FLOW'] = array([[1]], dtype=uint8)
data['SOLID'] = array([[2]], dtype=uint8)

image = imread(matlabPath + imageName + imageExtension)

maxY = len(image)
maxX = len(image[0])

material = zeros((maxX, maxY))
for y in range(maxY):
    line = ''
    for x in range(maxX):
        if str(type(image[y][x])) != "<class 'numpy.uint8'>":
            material[x][maxY - y - 1] = 2 if image[y][x][0] == 0 else 1
        else:
            material[x][maxY - y - 1] = 2 if image[y][x] == 0 else 1

        if (x % 10 == 0 and y % 10 == 0):
            line += '#' if material[x][maxY - y - 1] == 2 else ' '

    if (y % 10 == 0):
        print(line)

data['material'] = array(material)
savemat(matlabPath + imageName + '.mat', data)

print('Processing complete')
