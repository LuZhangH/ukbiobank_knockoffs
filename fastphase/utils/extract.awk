BEGIN {
  getline n < linesfile
  if(length(ERRNO)) {
    print "Unable to open linesfile '" linesfile "': " ERRNO > "/dev/stderr"
    exit
  }
}

NR == n { 
  print
  if(!(getline n < linesfile)) {
    if(length(ERRNO))
      print "Unable to open linesfile '" linesfile "': " ERRNO > "/dev/stderr"
    exit
  }
}
