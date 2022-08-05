import re

def naturalKey(decimals=False):

  if decimals:
    number_regex = "([\d]+(\.[\d]+)?)"
    numType = float
  else:
    number_regex = "([\d]+)"
    numType = int
    
  regex = re.compile(number_regex)
  
  def key(s):
    s = regex.split(s)
    
    l = []
    for i in s:
      try:
        l.append(numType(i))
      except (ValueError, TypeError):
        l.append(i)
    return l
  
  return key

integerKey = naturalKey(decimals=False)
decimalKey = naturalKey(decimals=True)
