import pickle

# read python dict back from the file
pkl_file = open('data/processed/labels/full_scores.pickle', 'rb')
mydict2 = pickle.load(pkl_file)
pkl_file.close()

print(mydict2)