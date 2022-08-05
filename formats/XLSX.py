import xlsxwriter
 
def exportCoords(filePath, posDict, coordsDict):
  
  wb = xlsxwriter.Workbook(filePath)

  for chromo in posDict:
    chromoPos = posDict[chromo]
    chromoCoords = coordsDict[chromo]
    nModels = len(chromoCoords)
    nCoords = len(chromoPos)
    
    ws = wb.add_worksheet()
    row = 0
    ws.write(row, 0, 'Chromosome:')
    ws.write(row, 1, chromo)
    
    row += 1
    ws.write(row, 0, 'Position')
    for j in range(nModels):
      col = j * 3
      j2 = j + 1
      ws.write(row, col+1, 'X %d' % j2)
      ws.write(row, col+2, 'Y %d' % j2)
      ws.write(row, col+3, 'Z %d' % j2)
    
    for j in range(nCoords):
      row += 1
      col = 0
      ws.write_number(row, col, chromoPos[j])
      
      for k in range(nModels):
        ws.write_number(row, col+1, chromoCoords[k,j,0])
        ws.write_number(row, col+2, chromoCoords[k,j,1])
        ws.write_number(row, col+3, chromoCoords[k,j,2])
        col += 3
  
  wb.close()
