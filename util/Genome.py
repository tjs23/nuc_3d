import os, sqlite3

try:
  # Python 3
  from urllib.request import urlopen
  from urllib.parse import urljoin
except ImportError:
  # Python 2
  from urllib2 import urlopen
  from urlparse import urljoin

from os.path import dirname, join

DB_FILE = join(dirname(dirname(__file__)), 'data', 'genome', 'genomeAssemblies.sqlite3')

def downloadFile(remoteFile, localFile, overwrite=False):
  
  if not overwrite and os.path.exists(localFile):
    return
  
  print('  Downloading %s to %s ...' % (remoteFile, localFile))
  
  response = urlopen(remoteFile)
  
  try:
    data = response.read()
    fileObj = open(localFile, 'wb')
    fileObj.write(data)
    fileObj.close()
    print(' ...done')
    
  except UnicodeDecodeError as err:
    print('  ...fail', err)

def getGenomeAssemblyReports(localDir='/data/assembly/',
                             url='ftp://ftp.ncbi.nlm.nih.gov',
                             dirPath='genomes/ASSEMBLY_REPORTS/All/',
                             fileExt='assembly.txt'):

  ftpPath = urljoin(url, dirPath)
  response = urlopen(ftpPath)
  data = response.read()
  
  print('Getting file list')
  
  filePaths = []
  for line in data.split('\n'):
    line = line.strip()
 
    if not line:
      continue
 
    fileName = line.split()[-1]
    
    if fileName.endswith(fileExt):
      filePaths.append(fileName, ftpPath + fileName)

  localPaths = []  
  for fileName, remotePath in filePaths:
    localPath = os.path.join(localDir, fileName)
    localPath = os.path.abspath(localPath)
    
    downloadFile(remotePath, localPath)
    localPaths.append(localPath)

  return localPaths
  
TABLE_SQL_1 = """
CREATE TABLE Assembly (
  id VARCHAR(24),
  name VARCHAR(12),
  organism TEXT,
  taxid INT,
  PRIMARY KEY (id)
);"""


TABLE_SQL_2 = """
CREATE TABLE Sequence (
  name VARCHAR(16),
  assembly VARCHAR(12) NOT NULL,
  seqName VARCHAR(16),
  gbName VARCHAR(16),
  seqLen INT,
  ucscName TEXT,
  FOREIGN KEY (assembly) REFERENCES Assembly(id),
  PRIMARY KEY (name)
);
"""

def getAssemblyData(assemblyFile, role='assembled-molecule'):

  fileObj = open(assemblyFile)
  
  name, organism, taxid = None, None, None
  
  chromos = []
  
  for line in fileObj:
    if line[0] == '#':
      if line.startswith('# Assembly Name:'):
        name = line.split()[-1]
        
      elif line.startswith('# Organism name:'):
        organism = ' '.join(line.split()[3:])

      elif line.startswith('# Taxid:'):
        taxid = line.split()[-1]
      
    else:
      data = line[:-1].split('\t')
      
      if len(data) > 1:
        if data[1] != role:
          continue
 
        try:
          seqName, seqRole, chromo, locn, genBank, rel, refSeq, asu, seqLen, ucscName = data
        except ValueError as err:
          print line
          raise(err)
        
        if refSeq == 'na':
          continue
        
        if seqLen == 'na':
          seqLen = 0
        
        if seqName == 'na':
          if ucscName != 'na':
            seqName = ucscName
          elif chromo != 'na':
            seqName = chromo + '?'
          else:
            seqName = '%d?' % (len(chromos)+1)
          
        chromos.append((refSeq, genBank, seqName, ucscName, int(seqLen)))
  
  fileObj.close()
  
  assembly = (name, organism, taxid)
  
  return assembly, chromos
  
HUMAN_NAMES = {'Drosophila_melanogaster': 'Fruit fly', 'Oryctolagus_cuniculus': 'Rabbit',
               'Fonticula_alba': 'Slime mould', 'Callorhinchus_milii': 'Shark','Nosema_bombycis': 'Microsporidian',
               'Monodelphis_domestica': 'Opossum', 'Volvox_carteri': 'Green alga', 'Sarcophilus_harrisii': 'Tasmanian devil',
               'Chelonia_mydas': 'Turtle', 'Bodo_saltans': 'Flagellate excavate', 'Giardia_lamblia': 'Parasitic excavate',
               'Moniliophthora_perniciosa': 'Rot fungus', 'Cryptococcus_neoformans': 'Dimorphic fungus',
               'Trypanosoma_brucei': 'Flagellate excavate', 'Gibberella_zeae': 'Blight fungus',
               'Dictyostelium_discoideum': 'Soil amoeba', 'Mortierella_verticillata': 'Soil fungus',
               'Pichia_sorbitophila': 'Yeast', 'Branchiostoma_floridae': 'Lancelet',
               'Coprinopsis_cinerea': 'Mushroom', 'Tuber_melanosporum': 'Truffle', 'Strigamia_maritima': 'Centipede',
               'Apis_mellifera': 'Honeybee', 'Latimeria_chalumnae': 'Coelacanth', 'Ornithorhynchus_anatinus': 'Platypus',
               'Phaeodactylum_tricornutum': 'Diatom', 'Trichomonas_vaginalis': 'Parasitic excavate',
               'Naegleria_gruberi': 'Excavate',  'Ostreococcus_tauri': 'Green alga',
               'Pythium_aphanidermatum': 'Oomycete', 'Crassostrea_gigas': 'Giant clam',
               'Monosiga_brevicollis': 'Choanoflagellate', 'Leishmania_major': 'Parasitic excavate', 'Physcomitrella_patens': 'Moss',
               'Anolis_carolinensis': 'Lizard', 'Ciona_savignyi': 'Tunicate', 'Batrachochytrium_dendrobatidis': 'Chytrid fungus',
               'Bigelowiella_natans': 'Rhizarian alga', 'Daphnia_pulex': 'Waterflea', 'Oryza_sativa': 'Rice',
               'Lottia_gigantea': 'Limpet', 'Salpingoeca_rosetta': 'Choanoflagellate', 'Dendroctonus_ponderosae': 'Beetle',
               'Allomyces_macrogynus': 'Chytrid fungus', 'Schizosaccharomyces_pombe': 'Yeast',
               'Toxoplasma_gondii': 'Parasitic alveolate', 'Takifugu_rubripes': 'Pufferfish',
               'Capitella_teleta': 'Polychaete', 'Tetrahymena_thermophila': 'Cilliate alveolate',
               'Rattus_norvegicus': 'Rat', 'Neurospora_crassa': 'Bread mould', 'Bombyx_mori': 'Silkmoth',
               'Strongylocentrotus_purpuratus': 'Sea urchin', 'Albugo_laibachii': 'Oomycete',
               'Rhizophagus_irregularis': 'Mycorrhizal fungus', 'Ustilago_maydis': 'Pathogenic fungus',
               'Sphaeroforma_arctica': 'Ichthyosporean', 'Hydra_vulgaris': 'Hydra', 'Perkinsus_marinus': 'Marine alveolate',
               'Puccinia_graminis': 'Rust fungus', 'Hyaloperonospora_arabidopsidis': 'Oomycete',
               'Plasmodium_falciparum': 'Parasitic alveolate',  'Acanthamoeba_castellanii': 'Amoeba',
               'Emiliania_huxleyi': 'Coccolithophore', 'Aspergillus_nidulans': 'Mould',
               'Colletotrichum_orbiculare': 'Pathogenic fungus', 'Canis_lupus': 'Dog', 'Clavispora_lusitaniae': 'Yeast',
               'Caenorhabditis_elegans': 'Nematode', 'Phytophthora_infestans': 'Oomycete',
               'Oncorhynchus_mykiss': 'Trout', 'Chlamydomonas_reinhardtii': 'Green alga', 'Capsaspora_owczarzaki': 'Opisthokont',
               'Xenopus_tropicalis': 'Frog', 'Felis_catus': 'Cat', 'Pristionchus_pacificus': 'Nematode',
               'Oikopleura_dioica': 'Tunicate', 'Petromyzon_marinus': 'Lamprey', 'Mus_musculus': 'Mouse', 'Saccharomyces_cerevisiae': 'Yeast',
               'Ectocarpus_siliculosus': 'Brown alga', 'Tetranychus_urticae': 'Spider mite', 'Paramecium_tetraurelia': 'Cilliate alveolate',
               'Sus_scrofa': 'Pig', 'Entamoeba_histolytica': 'Parasitic amoeboa', 'Trichoplax_adhaerens': 'Plate animal',
               'Arabidopsis_thaliana': 'Thale cress', 'Guillardia_theta': 'Cryptophyte alga', 'Reticulomyxa_filosa': 'Foraminifer',
               'Bos_taurus': 'Cow', 'Candida_albicans': 'Yeast', 'Homo_sapiens': 'Human',
               'Encephalitozoon_cuniculi': 'Parasitic fungus', 'Nematostella_vectensis': 'Sea anemone', 'Yarrowia_lipolytica': 'Yeast',
               'Thecamonas_trahens': 'Aspuzoan protist', 'Gallus_gallus': 'Chicken', 'Spizellomyces_punctatus': 'Chytrid fungus',
               'Anopheles_gambiae': 'Mosquito', 'Amphimedon_queenslandica': 'Sponge', 'Zea_mays': 'Maize',
               'Xiphophorus_maculatus': 'Platyfish', 'Thalassiosira_pseudonana': 'Marine diatom', 'Taeniopygia_guttata': 'Zebra finch',
               'Danio_rerio': 'Zebrafish', 'Selaginella_moellendorffii': 'Spikemoss', 'Loxodonta_africana':'Elephant'}

genera = set([x.split('_')[0] for x in HUMAN_NAMES.keys()])

def buildAssemblyDatabase(download=True, localDir='/data/assembly/', dbFile=DB_FILE):  
  localDir = os.path.abspath(localDir)
  
  if download:
    filePaths = getGenomeAssemblyReports(localDir)
  
  else:
    filePaths = []
    
    for fileName in os.listdir(localDir):
      localPath = os.path.join(localDir, fileName)
      filePaths.append(localPath)
  
  dbFile =  os.path.abspath(dbFile)
  
  if os.path.exists(dbFile):
    os.unlink(dbFile)
      
  connection = sqlite3.connect(dbFile)
  cursor = connection.cursor()
  cursor.execute(TABLE_SQL_1)
  cursor.execute(TABLE_SQL_2)
  cursor.execute('CREATE INDEX gb_index on Sequence (gbName)')
  cursor.execute('CREATE INDEX ucsc_index on Sequence (ucscName)')
  cursor.close()
  connection.commit()
 
  cursor = connection.cursor()
  done = set()
  
  for filePath in filePaths:
    assembly, chromos = getAssemblyData(filePath)
    
    if None in assembly:
      continue
    
    name, organism, taxid = assembly
    if organism.split()[0] not in genera:
      continue
       
    key = '%s_%s' % (name, taxid)
    
    #print name, taxid, organism, len(chromos)
    
    if key in done:
      print "* * * *  REPEAT  * * * *"
    
    else:    
      done.add(key)
      smt = 'INSERT INTO Assembly (id, name, organism, taxid) VALUES (?, ?, ?, ?)'
      
      try:
        cursor.execute(smt, (key, name, organism, taxid))
      except Exception  as err:
        cursor.close()
        connection.commit()
        print key, name, organism, taxid, type(key), type(name), type(organism), type(taxid)
        raise(err)
      
        
    for refSeq, genBank, seqName, seqLen, ucscName in chromos:
      if refSeq in done:
        print "* * * *  REPEAT  * * * *", refSeq
        continue
     
      done.add(refSeq)
      smt = 'INSERT INTO Sequence (name, assembly, seqName, gbName, ucscName, seqLen) VALUES (?, ?, ?, ?, ?, ?)'
      
      try:
        cursor.execute(smt, (refSeq, key, seqName, genBank, ucscName, seqLen))
      except Exception as err:
        print refSeq, filePath
        raise err
    
  cursor.close()
  connection.commit()
  connection.close()


def getChromosomeNames(chromoIds, dbFile=DB_FILE):

  dbFile =  os.path.abspath(dbFile)
  
  connection = sqlite3.connect(dbFile)
  cursor = connection.cursor()
  
  smt1 = 'SELECT seqName, ucscName FROM Sequence WHERE name=?'
  smt2 = 'SELECT seqName, ucscName FROM Sequence WHERE gbName=?'
  
  nameDict = {}
  
  for i, chromoId in enumerate(chromoIds):
    if chromoId[:3] == 'gi|':
      query = chromoId.split('|')[3]
    else:
      query = chromoId
        
    name = chromoId
    
    for smt in (smt1, smt2):
      result = cursor.execute(smt, [query]).fetchone()
      
      if result:
        name, ucscName = result
        break
    
    nameDict[chromoId] = name
    
  cursor.close()
  connection.close()

  return nameDict
  

if __name__ == '__main__':

  buildAssemblyDatabase()













