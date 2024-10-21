
def check_eligibility(counts,MinClusterSize):
    if max(counts)> sum(counts)/MinClusterSize:
        return False
    return True


def alpha_beta(df,sensitive_column,MinClusterSize):
    #compute alpha and beta for remainig tuples without cluster
    value_counts = df.loc[df['cluster'] == 0, sensitive_column].value_counts().to_list()
    remaining = sum(value_counts)
    alpha = 1
    beta = MinClusterSize

    #find first alpha that works
    found = False
    while(not found):
        tmp = (remaining-alpha*beta)/MinClusterSize
        if beta < len(value_counts):
            if(value_counts[0]-alpha <= tmp and value_counts[beta]<= tmp and alpha<=value_counts[beta-1]):
                found = True
            else:
                beta +=1
        elif beta == len(value_counts):
            if(value_counts[0]-alpha <= tmp and alpha<=value_counts[beta-1]): 
                found = True
        else:
            print("Something went wrong, beta> num values")
    
    #find biggest alpha that works
    tmp = (remaining-alpha*beta)/MinClusterSize
    
    if(beta<len(value_counts)):
        while((value_counts[0]-alpha) <= tmp and (value_counts[beta]<= tmp) and (alpha<=value_counts[beta-1])):
            alpha += 1
            tmp = (remaining-alpha*beta)/MinClusterSize
    else:
        while( value_counts[0]-alpha <= tmp and alpha<=value_counts[beta-1]):
            alpha += 1
            tmp = (remaining-alpha*beta)/MinClusterSize
    
    #fix alpha to last functional value
    alpha -= 1
    return(alpha,beta)







          
           



def compute_clusters(df,sensitive_column,MinClusterSize):
    #Assigns tuples to cluster on column "clusters"
    

    cluster_number = 1

    value_counts = df.loc[df['cluster'] == 0, sensitive_column].value_counts().to_list()
    value_names  = df.loc[df['cluster'] == 0, sensitive_column].value_counts().index.to_list()
    #print(value_counts)

    while(len(value_names)>1): #while available tuples


        if not check_eligibility(value_counts,MinClusterSize):
            print("Not m-eligible")
            print(value_counts)
            return(df)
        alpha,beta = alpha_beta(df,sensitive_column,MinClusterSize)
        #print(alpha,beta)

        #With alpha and beta computed we create clusters. They are made greedily

        #create alpha clusters with size beta

        #for each attribute 
        for attribute in value_names[:beta]:
            selected_rows = df.loc[(df['cluster'] == 0) & (df[sensitive_column]==attribute)]

            if not len(selected_rows)<alpha:
                first_alpha_occurrence_index = selected_rows.index[:alpha]
                df.loc[first_alpha_occurrence_index, 'cluster'] = list(range(cluster_number,cluster_number+alpha))
            else:
                print('Error: len(selected_rows)<alpha')
        cluster_number += alpha

        value_counts = df.loc[df['cluster'] == 0, sensitive_column].value_counts().to_list()
        value_names  = df.loc[df['cluster'] == 0, sensitive_column].value_counts().index.to_list()
        #print(value_counts)
        #print(value_names)


        
    if min(df["cluster"])==0:
        print("There are tuples with non assigned cluster!")
    return(df)





def optimize_clusters(df,sensitive_column, quasiidentifiers,generalize=False):
    import numpy as np
    import scipy.spatial.distance as d
    #greedy method to make better clusters. TO DO
    #search all clusters with common signature, then regroup tuples to minimize some value function
    signature_dict = {}
    cluster_numbers  = df["cluster"].value_counts().index.to_list() #probably equal to max but just in case
    #(cluster_numbers)
    
    #merge clusters with common signature
    for cl_number in cluster_numbers:
        #merge clusters with common signature
        signature_cl = str(sorted(list(df.loc[df['cluster'] == cl_number, sensitive_column].unique()))) #store signature as a string, sort to make strings equal
        
        if signature_cl in signature_dict.keys():
            df.loc[df['cluster'] == cl_number, "cluster"] = signature_dict[signature_cl]
        else:
            signature_dict[signature_cl] = cl_number #use same number to avoid repetition of values with other clusters
    
    #return(signature_dict)
    
    #split clusters with common signature heuristic
    cluster_numbers  = df["cluster"].value_counts().index.to_list() #update cluster identifiers
    #good step to divide shared identifiers...

    
    #for each bucket must be done alpha times!
    new_cl_index =-1

    for cl_number in cluster_numbers:
        #while tuples in bucket

        while(len(df.loc[df["cluster"]==cl_number])>0):
            #with signature...
            signature_list = list(df.loc[df['cluster'] == cl_number, sensitive_column].unique())
            #has tuples...
            #selected_rows = df.loc[df['cluster'] == cl_number] #select spefic tuples to organize

            #pick an initial tuple with sensitive value signature_list[0]
            first_tuple_index = df.loc[(df['cluster'] == cl_number) & (df[sensitive_column]==signature_list[0])].index[0]
            df.loc[first_tuple_index, 'cluster'] = new_cl_index
            #compute centroid of cluster -> asume numerical values from now on
            #centroid = average QI of tuples in cluster
            centroid = df.loc[df["cluster"]== new_cl_index,quasiidentifiers].mean() ##!!! have to remove sensitive columns or specify QI!
            #print(centroid)
            best_dist = np.inf
            for attribute in signature_list[1:]:
                #pick closest to cluster to add with attribute
                candidates = df.loc[(df['cluster'] == cl_number) & (df[sensitive_column] == attribute)]  # & (selected_rows['cluster']>0 )] #use only tuples of 
                for candidate in candidates.iterrows():
                    #return(centroid,candidate)
                    #print(centroid)
                    #print(candidate[1][quasiidentifiers])
                    q_values = np.array(candidate[1][quasiidentifiers],dtype=float)
                    if d.pdist([centroid,q_values],'euclidean')<best_dist:
                        best_dist = d.pdist([centroid,q_values],'euclidean')
                        best_candidate_index = candidate[0]
                
                best_dist = np.inf
                df.loc[best_candidate_index,"cluster"] = new_cl_index
                centroid = df.loc[df["cluster"]== new_cl_index,quasiidentifiers].mean()
                
                
                #update cluster with chosen candidate
                #centroid #HERE
            new_cl_index-= 1
    df['cluster'] = -df["cluster"]
    if(generalize):
        for x in df["cluster"].value_counts().index.to_list():
            df.loc[df["cluster"]==x,quasiidentifiers] = list(df.loc[df["cluster"]==x,quasiidentifiers].mean())
    return(df)

def randomize_clusters(df,sensitive_column,method='choice'):
    import random
    #set random seed 
    random.seed(2024)
    #NumClusters = max(df["cluster"])
    cluster_numbers  = df["cluster"].value_counts().index.to_list() #probably equal to max but just in case

    for x in cluster_numbers:
        value_names  = df.loc[df['cluster'] == x, sensitive_column].value_counts().index.to_list()
        #compute overwriting variation and write it on top of 
        if method == 'choice':
            variation = random.choices(value_names, k=len(value_names)-1) #with replacement -> at least 1 replacement
            repeated_element = random.choice(variation)

            variation_with_repetition = variation + [repeated_element]
            random.shuffle(variation_with_repetition) 
        elif method == "sample":
            variation = random.sample(value_names, k=len(value_names)-1) #without replacement -> forces exactly 1 replacement
            repeated_element = random.choice(variation)

            variation_with_repetition = variation + [repeated_element]
            random.shuffle(variation_with_repetition) 
        elif method == "ord_sample":
            #sample but keeping the order
            random_index = random.randint(0, len(value_names) - 1)
            random_index2 = random.randint(0, len(value_names) - 2)
            variation_with_repetition = value_names.copy()
            if random_index2<random_index:
                variation_with_repetition[random_index] = variation_with_repetition[random_index2]
            else:
                variation_with_repetition[random_index] = variation_with_repetition[random_index2+1]

            #variation = random.sample(value_names, k=len(value_names)-1) #without replacement -> forces exactly 1 replacement
            #repeated_element = random.choice(variation)

            #variation_with_repetition = variation + [repeated_element]
            #random.shuffle(variation_with_repetition) 
        elif method == "permutation": #simply permute the values of each cluster
            variation_with_repetition = random.sample(value_names, k=len(value_names)) # permutation -> no replacement
            random.shuffle(variation_with_repetition) 
        else:
            print("Invalid method")

        df.loc[df['cluster'] == x, sensitive_column] = variation_with_repetition

    return(df)

from scipy.spatial.distance import hamming

def dist(tuple1,tuple2):
    return hamming(tuple1,tuple2) #other alternatives are possible

def total_distance(df):
    total_dist=0
    unique_clusters = df["cluster"].unique()

    for cluster in unique_clusters:
        cluster_tuples = df.loc[df["cluster"]==cluster]
        size = len(cluster_tuples)
        total_dist += sum(dist(cluster_tuples.iloc[i],cluster_tuples.iloc[j]) for i in range(size) for j in range(i+1,size))
    
    return(total_dist)

def check_swap(df,index1,index2):
    cluster_1,cluster_2 = df.iloc[index1]["cluster"],df.iloc[index2]["cluster"]
    old_dist = 0
    for cluster in [cluster_1,cluster_2]:
        cluster_tuples = df.loc[df["cluster"]==cluster]
        size = len(cluster_tuples)
        old_dist += sum(dist(cluster_tuples.iloc[i],cluster_tuples.iloc[j]) for i in range(size) for j in range(i+1,size))
    
    new_dist = 0

    cluster_tuples = df.loc[(df["cluster"]==cluster_1) & (df.index != index1)]
    size = len(cluster_tuples)
    new_dist += sum(dist(df.iloc[index2],cluster_tuples.iloc[j]) for j in range(size)) #all tuples of cluster with new one
    new_dist += sum(dist(cluster_tuples.iloc[i],cluster_tuples.iloc[j]) for i in range(size) for j in range(i+1,size)) #tuples of cluster with each other

    cluster_tuples = df.loc[(df["cluster"]==cluster_2) & (df.index != index2)]
    size = len(cluster_tuples)
    new_dist += sum(dist(df.iloc[index1],cluster_tuples.iloc[j]) for j in range(size)) #all tuples of cluster with new one
    new_dist += sum(dist(cluster_tuples.iloc[i],cluster_tuples.iloc[j]) for i in range(size) for j in range(i+1,size)) #tuples of cluster with each other

    return(new_dist<old_dist,new_dist,old_dist)


def do_swap(df,index1,index2):
    df.loc[df.index == index1,"cluster"],df.loc[df.index == index2,"cluster"] = df.loc[df.index == index2,"cluster"],df.loc[df.index == index1,"cluster"]
    return


from numpy import Inf

def greedy_swap(df,shared_identifiers,sensitive_column):
    improvement_val = 0
    best_dist = Inf
    improvement =True
    #while(improvement):
    improvement= False
    for index1 in range(1000): #range(len(df)):
        improvement = False
        #cluster of index1
        cluster1 = df.iloc[index1]["cluster"]
        #is index1 a shared_identifier
        shared1 = df.iloc[index1]["identifier"] in shared_identifiers
        #compute number of shared identifiers in cluster
        #shared_c1 = len(df.loc[df["cluster"] == cluster1 & (df["identifier"].isin(shared_identifiers))])
        signature1 = df.loc[df["cluster"]==cluster1,sensitive_column].to_list()
        for index2 in range(index1+1,1000):  #len(df)):
            #cluster of index2
            cluster2 = df.iloc[index2]["cluster"]
            if cluster1 != cluster2:
                #is index2 a shared_identifier
                shared2 = df.iloc[index2]["identifier"] in shared_identifiers
                
                if df.iloc[index1][sensitive_column] == df.iloc[index2][sensitive_column]: #same attribute implies different cluster
                    swap,new_dist,old_dist = check_swap(df,index1,index2)
                    if swap:
                        #compute number of shared identifiers in cluster
                        #shared_c2 = len(df.loc[df["cluster"] == df.iloc[index2]["cluster"] & (df["identifier"].isin(shared_identifiers))])
                        #condition = max(shared_c1-shared1+shared2,shared_c2-shared2+shared1)<max(shared_c1,shared_c2)
                        condition = shared1==shared2
 
                        if condition and new_dist <best_dist: #check if it causes problems with shared identifiers, if fine, do swap
                            #do_swap(df,index1,index2)
                            best_dist = new_dist
                            best_old_dist = old_dist
                            best_index = index2
                            improvement =True

                

                else:  #swapping preserves m-uniqueness
                    signature2 = df.loc[df["cluster"]==cluster2,sensitive_column].to_list()
                    if df.iloc[index1][sensitive_column] not in signature2 and df.iloc[index2][sensitive_column] not in signature1:
                        swap,new_dist,old_dist = check_swap(df,index1,index2)
                        if swap:
                            #compute number of shared identifiers in cluster
                            #shared_c2 = len(df.loc[df["cluster"] == df.iloc[index2]["cluster"] & (df["identifier"].isin(shared_identifiers))])
                            condition = shared1 == shared2 #max(shared_c1-shared1+shared2,shared_c2-shared2+shared1)<max(shared_c1,shared_c2)
                            if condition and new_dist <best_dist: #check if it causes problems with shared identifiers, if fine, do swap
                                #do_swap(df,index1,index2)
                                best_dist = new_dist
                                best_old_dist = old_dist
                                best_index = index2
                                improvement =True
        if(improvement):
            do_swap(df,index1,best_index)
            improvement_val = improvement_val+best_dist-best_old_dist
    return(improvement_val)

                
                    






