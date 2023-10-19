from Bio import Align
from Bio.Seq import Seq
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from sklearn.cluster import KMeans
from sklearn.metrics import accuracy_score
import random

aligner = Align.PairwiseAligner()
dataset = pd.read_csv('dataset.txt',nrows=100)
dataset['querydate']= pd.to_datetime(dataset['querydate'])



train = dataset[['latitude','longitude']]
kmeans = KMeans(n_clusters=2)
kmeans.fit(train)
predict = kmeans.predict(train)


cluster_labels = kmeans.labels_
for i in range(0,10):
    predict[i] = 3
acc = accuracy_score(cluster_labels,predict)
print(acc)
dataset['clusterID'] = cluster_labels
print(dataset.head())

dataset = dataset.values
trajectory_append = []
loss = []
hit = 0

def dynamicSA(src_lat,src_lon,cls_id):
    global hit
    dups = []
    max1 = 0
    max2 = 0
    choosen_lat = 0
    choosen_lon = 0
    while len(dups) < len(dataset):
        random_record_dataset = 0
        flag = True
        while flag:
            random_record_dataset = random.randint(0,(len(dataset)-1))
            if random_record_dataset not in dups:
                dups.append(random_record_dataset)
                flag = False            
        des_lat = dataset[random_record_dataset,2]
        des_lon = dataset[random_record_dataset,3]
        seq1 = Seq(str(src_lat))
        seq2 = Seq(str(des_lat))
        seq3 = Seq(str(src_lon))
        seq4 = Seq(str(des_lon))
        alignments1 = aligner.align(seq1, seq2)
        alignments2 = aligner.align(seq3, seq4)
        if not alignments1:
            continue
        for match in alignments1:
            score = match.score
            if score > max1:
                max1 = score
                choosen_lat = des_lat
        if not alignments2:
            continue        
        for match in alignments2:
            score = match.score
            if score > max2:
                max2 = score
                choosen_lon = des_lon
    cls = 0
    if max1 <= 5 and max2 <= 5:
        cls = 0
    else:
        cls = 1
    if cls == cls_id:
        hit = hit + 1
    print(str(hit)+" "+str(cls_id)+" "+str(max1)+" "+str(max2))    
    return str(choosen_lat)+","+str(choosen_lon),score  
    
total = 0
for i in range(len(cluster_labels)):
    j = i + 1
    max1 = 0
    max2 = 0
    choosen_lat = 0
    choosen_lon = 0
    src_lat = dataset[i,2]
    src_lon = dataset[i,3]
    cls_id = dataset[i,4]
    trajectory_value, trajectory_loss = dynamicSA(src_lat,src_lon,cls_id)
    trajectory_append.append(trajectory_value)
    loss.append(trajectory_loss)
    print(trajectory_value+" "+str(trajectory_loss))
    
print(hit)   