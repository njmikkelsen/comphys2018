import sys, os
import requests

def download_raw_file(url,filename):
  r = requests.get(url)
  with open(filename, 'w') as F:
    F.write(r.text)


def verify_files():
  urlbase = "https://raw.githubusercontent.com/njmikkelsen/comphys2018/master/Project1/"
  repo    = "https://github.com/njmikkelsen/comphys2018/tree/master/Project1"
  print("Verifying files...")
  while True:
    try:
      filename = "tridiagmateq.cpp"; open(filename)
      filename = "tridiagmateq.h";   open(filename)
      break
    except:
      print("Error: Cannot find file '{:s}' in the local directory!".format(filename))
      print("Do you want to download the file? (yes/no)"); k = 0
      while True:
        user_input = input(" | ")
        if user_input.lower() == "yes":
          print("Downloading file...")
          download_raw_file(urlbase+filename,filename)
          print("Done.")
          break
        elif user_input.lower() == "no" or k > 1:
          if k > 1: print("Too many attempts.")
          print("You may manually download the file from the following repository:")
          print(repo)
          print("Unable to continue without the file, quitting the program.")
          sys.exit(1)
        else: k += 1
  if not os.path.isdir("input"):
    print("Missing directory 'input'")
    os.mkdir("input")
  if not os.path.isdir("output"):
    print("Created directory 'output'")
    os.mkdir("output")
  print("All files accounted for.")











