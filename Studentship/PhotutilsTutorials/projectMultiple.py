from astropy.io import fits
from reproject import reproject_interp, reproject_adaptive, reproject_exact
from astropy.wcs import WCS
import matplotlib.pyplot as plt
#objectName = 'JL163'
#filterType = 'V'
#originalFileName = 'Studentship/SPIRIT6-Observations/GoodDarks/Multiple-Exposures/JL-163/JL-163-S001-V-R001-V.fts'

def imageProject(objectName, filterType, originalFileName, date:str, mode):
    """
    if mode == 'single': \n
        objectName - the name of an individual object, e.g. 'JL163' \n
        filterType - letter corresponding to filter, 'B' or 'V' \n
        originalFileName - this is the base image that all other images will be projected onto. give the path to the image. \n
        date - yyyymmdd format e.g. '20230112' for Jan 12th 2023. \n
    \n 
    elif mode == 'multiple': \n
        date - yyyymmdd,yyyymmdd.... repeated for as many datasets as you want of the given object e.g. 20230112,20230114,20230118... \n
    """
    if mode == 'single':

        file = f'Studentship/PhotutilsTutorials/90secExposures{objectName}_{filterType}_{date}.txt'
        f2 = open(f'Studentship/PhotutilsTutorials/{date}_{objectName}ReprojectedImages-{filterType}.txt', 'w')
        f2.write(originalFileName + '\n')
        f2.close()
        with open(file, 'r') as f:
            lines = f.readlines()
            init = lines[0].strip('\n')
            hdu1 = fits.open(init)[0]     
            for line in lines[1:]:
                line2 = line.strip('\n')
                var = line.strip('\n').split('-')
                if objectName == 'JL163':
                    num = var[7]
                elif objectName == 'WD0830':
                    num = var[9]
                elif objectName == 'LB1735':
                    num = var[7]
                elif objectName == 'MCT0401':
                    num = var[9]
                elif objectName == 'MCT0550':
                    num = var[9]
                elif objectName == 'TPhe':
                    num = var[7]
                elif objectName == 'WD1153':
                    num = var[9]
                hdu2 = fits.open(line2)[0]
                array, footprint = reproject_exact(hdu2, hdu1.header)
                fits.writeto(f'Studentship/ReprojectedImages/{date}_{objectName}Reprojected-{num}-{filterType}.fits', array, hdu1.header, overwrite=True)
                f3 = open(f'Studentship/PhotutilsTutorials/{date}_{objectName}ReprojectedImages-{filterType}.txt', 'a')
                f3.write(f'Studentship/ReprojectedImages/{date}_{objectName}Reprojected-{num}-{filterType}.fits\n')
                f3.close()
    if mode == 'multiple':
        dateList = date.split(',')
        originalDate = dateList[0]

        hdu1 = fits.open(originalFileName)[0]
        hdu2List = []

        for item in dateList:
            file = f'Studentship/PhotutilsTutorials/90secExposures{objectName}_{filterType}_{item}.txt'
            with open(file, 'r') as f:
                lines = f.readlines()    
                for line in lines:
                    line2 = line.strip('\n')
                    var = line.strip('\n').split('-')
                    hdu2 = fits.open(line2)[0]
                    hdu2List.append(hdu2)


        f3 = open(f'Studentship/PhotutilsTutorials/{date}_{objectName}ReprojectedImages-{filterType}.txt', 'a')   
        num = 1         
        for hdu in hdu2List:
            array, footprint = reproject_exact(hdu, hdu1.header)
            fits.writeto(f'Studentship/ReprojectedImages/{date}_{objectName}Reprojected-{num}-{filterType}.fits', array, hdu1.header, overwrite=True)
            f3.write(f'Studentship/ReprojectedImages/{date}_{objectName}Reprojected-{num}-{filterType}.fits\n')
            num += 1
        f3.close()

#imageProject('WD1153', 'V', 'Studentship/SPIRIT6-Observations/GoodDarks/Multiple-Exposures/20230118/WD-1153-484/WD-1153-484-S001-V-R001-V.fts', '20230112,20230114,20230118','multiple')
        




