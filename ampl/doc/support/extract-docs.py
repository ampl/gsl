#!/usr/bin/env python
"""Extract documentation from source file.
Usage: extract-docs.py <sourcefile> [-o DIR  Output directory for generated rst files]
"""
import os, re
from docopt import docopt

def extract_docs(filename, output_dir):
  "Extract the AMPLGSL documentation from the code."
  output = None
  if not os.path.exists(output_dir):
    os.mkdir(output_dir)
  with open(filename, 'r+t') as input:
    map = input.read()
    for i in re.finditer(r'/\*\*(.*?)\*/', map, re.DOTALL):
      s = re.sub(r'\n +\* ?', r'\n', i.group(1))
      s = re.sub(r'\$(.+?)\$', r':math:`\1`', s, flags=re.DOTALL)
      m = re.search(r'@file (.*)', s)
      if m:
        filename = m.group(1)
        if output:
          output.close()
        output = open(os.path.join(output_dir, filename + '.rst'), 'w')
        s = s[:m.start()] + s[m.end():]
      output.write(s.rstrip(' '))

if __name__ == '__main__':
  args = docopt(__doc__)
  if(args['-o']):
    outputdir = args['DIR']
  else:
    outputdir = './amplgsl'
  print(f"Extracting {args['<sourcefile>']} to {outputdir}")
  extract_docs(args['<sourcefile>'], outputdir)
