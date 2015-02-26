import xml.etree.ElementTree as ET
import re
from collections import Counter
import os
import json
import sys
from math import sqrt, log

from stemming.porter2 import stem

from pprint import pprint
from flask import Flask, redirect, request

app = Flask(__name__)



root_dir = "E:/LocalData/PubMed/Cell/"

stop_words = [w.replace("\n", "") for w in open("stop_words.txt")]

unstemmed = {}

def file_root(f):
    tree = ET.parse(f)
    return tree.getroot()

def get_abstract_node(root):
    return root.iter("abstract").next()

def get_text(node):

    ts = [get_text(n) for n in node]
    
    txts = [node_text for sublist in ts for node_text in sublist]

    if node.text:
        txts.append(node.text)
    return txts

def is_float(f):
    try:
        float(f)
        return True
    except:
        return False

def words(text_list):
    word_list = [re.split('\W+',t) for t in text_list]
    
    words = [unicode(word).lower() for text in word_list for word in text]
    
    # Filter out words that are too short.
    words = [w for w in words if len(w) > 2]
    
    # Filter out words that are just numbers.
    words = [w for w in words if not is_float(w)]
    
    # Filter out words that start with a number.
    words = [w for w in words if not is_float(w[0])]
    
    # Filter out stop words
    words = [w for w in words if not w in stop_words]
    
    # Do stemming
    stemmed = [(w, stem(w)) for w in words]
    
    # Save map of stem -> original
    for w,s in stemmed:
        unstemmed[s] = w
    
    return [s for w,s in stemmed]

def normalise_counter(c):
    total_words = sum(c.values())
    
    for a in c:
        c[a] /= float(total_words)

def process_file(f):
    
    root = file_root(f)
    
    try:
        abstract = get_abstract_node(root)
    except:
        return Counter()
    
    all_words = words(get_text(abstract))
    
    return Counter(all_words)

def process_all_files():
    
    corpus_words = Counter()
    
    files = [os.path.join(root_dir, f) for f in os.listdir(root_dir)]

    word_files = {}
    file_words = {}
    for f in files:
        print "Processing %s" % f
        sys.stdout.flush()
        
        file_words[f] = process_file(f)
        corpus_words.update(file_words[f])
        for k in file_words[f]:
            word_files[k] = word_files.get(k, set())
            word_files[k].add(f)
            
    return file_words, corpus_words, word_files, files    


file_words, corpus_words, word_files, files = process_all_files()

print "Computed word count for entire corpus"




## K-Means stuff

def distance(fv1, fv2):
    all_features = set([f for f in fv1]).union([f for f in fv2])
    
    s = 0
    for f in all_features:
        s += (fv1.get(f,0) - fv2.get(f,0)) ** 2
        
    return sqrt(s)

def centroid(fvs):
    
    centre = Counter()
    
    for f in fvs:
        centre += f
        
    for a in centre:
        centre[a] /= float(len(fvs))
        
    return centre

def nearest_centroid_index(fv, centroids):
    
    min_dist = distance(fv, centroids[0])
    nearest_index = 0
    
    for i in range(1, len(centroids)):
        d = distance(fv, centroids[i])
        if d < min_dist:
            min_dist = d
            nearest_index = i
    
    return nearest_index, min_dist

def kmeans(xs,K):
    """xs is a vector of Counters. K is the required number of clusters"""

    print "DISTANCE: " + str(distance(xs[0], xs[1]))
    print xs[0]
    print xs[1]
    print xs[2]

    N = len(xs)
    seeds = xs[:K] # Pick first K papers as seed cluster centroids

    # Initialise centroids
    mu = seeds

    # Initialise cluster members
    omega = []
    for i in range(K):
        omega.append(set())


    last_omega = omega[:]
    i = 0
    while i < 20:
        for k in range(K):
            omega[k] = set()

        # Reassign xs to nearest centroid

        for n in range(N):
            j,dist = nearest_centroid_index(xs[n], mu)
            omega[j].add(n)
            print "Assigning paper %d to cluster %d (distance %f)" % (n,j,dist)

        print "Omega " + str(i) + ": " + str(omega)

        # Recompute cluster centres

        for k in range(K):
            mu[k] = centroid([xs[t] for t in omega[k]])

        match = True
        for t in range(len(omega)):
            if omega[t] != last_omega[t]:
                match = False
                break

        if match:
            print "\nConverged in %d iterations.\n\n" % i
            break 


        i += 1

        last_omega = omega[:]
    
    return omega
    

def mutual_information(word, cluster, files, feature_vectors):
    N = float(len(files))
    
    # The number of files not containing the word and not in the cluster
    N00 = float(len([f for fi,f in enumerate(files) if not fi in cluster and not word in feature_vectors[fi]]))
    
    # The number of files containing the word and in the cluster
    N11 = float(len([f for fi,f in enumerate(files) if fi in cluster and word in feature_vectors[fi]]))
    
    # The number of files containing the word, but not in the cluster
    N10 = float(len([f for fi,f in enumerate(files) if not fi in cluster and word in feature_vectors[fi]]))
    
    # The number of files not containing the word, but in the cluster
    N01 = float(len([f for fi,f in enumerate(files) if fi in cluster and not word in feature_vectors[fi]]))
    
    # The number of files containing the word
    N1_ = N11 + N10
    
    # The number of files not containing the word
    N0_ = N00 + N01
    
    # The number of files in the cluster
    N_1 = N01 + N11
    
    # The number of files not in the cluster
    N_0 = N00 + N10
    
    #print N00, N01, N10, N11
    
    mi = 0
    if N11 > 0:
        mi += N11/N * log(N*N11/(N1_*N_1),2)
    if N01 > 0:
        mi += N01/N * log(N*N01/(N0_*N_1),2)
    if N10 > 0:
        mi += N10/N * log(N*N10/(N1_*N_0),2)
    if N00 > 0:
        mi += N00/N * log(N*N00/(N0_*N_0),2)
    
    return mi
         
    
def cluster_label(cluster, feature_vectors, files):

    if len(cluster) < 1:
        return "EMPTY"

    cluster_words = Counter()
    for fi in cluster:
        cluster_words += feature_vectors[fi]
        
    for w in cluster_words:
        cluster_words[w] = mutual_information(w, cluster, files, feature_vectors)
    
    return cluster_words.most_common(1)[0][0]
            
          


@app.route("/")
def root():
	return redirect("/static/index.html")


@app.route("/kws")
def kws():

	try:
		filter_kws = request.args['filter'].split(",")
	except:
		filter_kws = []

	# Only allow words that occur in at least two papers

	allowed_words = [w for w in word_files if len(word_files[w]) >= 2 and not w in filter_kws]

	filtered_files = []
	for f in files:
		include_file = len(filter_kws) == 0
		for w in file_words[f]:
			if w in filter_kws:
				include_file = True

		if include_file:
			filtered_files.append(f)

	print "FILTERED FILES:"
	print len(filtered_files)

	feature_vectors = []
	for f in filtered_files:
	    counter = Counter()
	    for w in file_words[f]:
	        if w in allowed_words:
	            counter[w] = file_words[f][w]
	            
	    feature_vectors.append(counter)


	clusters = kmeans(feature_vectors, min(50, len(filtered_files)))

	labels = []
	for cluster in clusters:
	    #if len(cluster) < 2:
	    #    continue # Ignore singleton clusters
	        
	    label = cluster_label(cluster, feature_vectors, filtered_files)
	    labels.append(label)
	    
	#pprint(labels)

	def files_with_both_keywords(kw1, kw2):
	    return word_files[kw1].intersection(word_files[kw2])


	def calculate_kw_matrix():

	    kw_matrix = {}
	    kw_papers = {}

	    for kw1 in all_kws:
	        print "Testing %s" % kw1
	        sys.stdout.flush()

	        kw_papers[kw1] = {
	            "papers": len(word_files[kw1])
	        }
	        kw_matrix[kw1] = {}

	        for kw2 in all_kws:
	            if kw2 != kw1:
	                papers_with_both = len(files_with_both_keywords(kw1,kw2))
	                if papers_with_both >= 1:
	                	kw_matrix[kw1][kw2] = papers_with_both

	        #if len(kw_matrix[kw1]) == 0:
	            #del kw_matrix[kw1]
	            #del kw_papers[kw1]
	            
	    return kw_matrix, kw_papers

	all_kws = labels

	kw_matrix, kw_papers = calculate_kw_matrix()

	print "Done. There are %d keywords in matrix." % len(kw_matrix)	   

	r = {
		"matrix": kw_matrix,
		"papers": kw_papers
	}
	return json.dumps(r, sort_keys=True, indent=2, separators=(',', ': '))

app.run(debug=True)