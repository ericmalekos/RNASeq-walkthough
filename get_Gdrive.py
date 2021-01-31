'''
This script serves as a template for copying multiple files from Google Drive
to the location the script is being run from.

** IMPORTANT ** Set the permission of the file in GDrive to "Anyone with link"

Dependencies:
	gdown

install with: pip3 install gdown --user

Variables:
	url_prefix (str): prefix of the GDrive share link should be the same in all cases
	suffix (str): the file type being transferred. Here we are assuming all file types
					are ".fg.gz".
	filedict (dict): The key is the name to save the file as, the value is the
					important portion of the link from GDrive.
					EX:
					Link: https://drive.google.com/file/d/1gc0Nfl693O49BECq34g6zno-ThJdBqw_/view?usp=sharing
					Relevant Portion: 1gc0Nfl693O49BECq34g6zno-ThJdBqw_
'''
import gdown
url_prefix = "https://drive.google.com/uc?id="
suffix = ".fq.gz"

filedict = {"A01_1" : "1DBFdAhsajwn25aUDdTJZbkANSfs63",
			"A01_2" : "1BdQWtJbk1Vj4EtGriXBGGHIBkPYsh",
			"A02_1" : "1el_4r7NdntMPP9qCPYk3s_YvWvZgx",
   			"A02_2" : "1lpwEzkV3y29SzLyGJB64tozhvs0V_",
			"A03_1" : "1uq8eMklE2v5YjtgxD1DutFmoB5Zcz",
   			"A03_2" : "1V7-bt2ECVrAQnbK3UO59JP1sJRlmz",
			"A04_1" : "1KFu7k6uq_InXbgjk-E6COUCWTlNY9",
   			"A04_2" : "1rtWLCAKzTsb-9hVU90bzHnc_n0Q0A",
			"A05_1" : "1TQ6DbO0-ki9vA-X-YlE2w5uCCHpSf",
			"A05_2" : "1_GdPa7RpilWLml8Z3Zgj-KKlIXP-X"}

for key, value in filedict.items():
	gdown.download(url_prefix + value, key + suffix, quiet=False)
