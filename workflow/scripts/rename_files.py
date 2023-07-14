#!/usr/bin/env python


import os
import re



aFiles = os.listdir()

for szFile in aFiles:
    # old name looks like:
    # stdin.part_093.fastq.gz
    # and we want to change it to:
    # 93-of-100.fq.gz

    szTemp   = re.sub( r'stdin.part_', '', szFile )
    szNumber = re.sub( r'.fastq.gz', '', szTemp )

    nNumber = int( szNumber )

    szNewName = str( nNumber ) + "-of-100.fq.gz"

    print( "renaming " + szFile  + " to " + szNewName )
    os.rename( szFile, szNewName )

    
    


