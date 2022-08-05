import sys

PYTHON_MAJOR = sys.version_info[0]

"""
`str` is for text data only, therefore `str` in py3 and `unicode` in py2.
`bytes` is for binary data only.
`tobytes` helper defines a consistent method to get a `bytes` object that will handle `str` inputs.
"""

if PYTHON_MAJOR == 3:
  str = str
  tobytes = lambda *args, **kwargs: bytes(*args, **kwargs)

elif PYTHON_MAJOR == 2:
  str = unicode

  def tobytes(x, encoding=None):
    if isinstance(x, (str, unicode)):
      if not encoding:
        raise TypeError("String type without encoding")
      return x.encode(encoding)
    else:
      # bytes is just an alias for `str` in py2
      return bytes(x)

else:
  raise RuntimeError("Unsupported Python version (%i)." % PYTHON_MAJOR)
