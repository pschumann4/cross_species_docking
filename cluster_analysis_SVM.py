import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.cluster import KMeans
import seaborn as sns
from kneed import KneeLocator
from sklearn import metrics
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
from sklearn import svm
from sklearn.model_selection import train_test_split

def cluster_analysis():
    
    # Read the summary csv file
    summary_file = input("Enter the file path for the summary file: ")
    # Remove the quotes from the file path
    summary_file = summary_file.replace('"', "")
    file_dir = os.path.dirname(summary_file)
    df = pd.read_csv(summary_file)
    # Multiply the "binding_affinity" and "lig_rmsd" columns by -1 so that higher values indicate a better pose
    df['binding_affinity'] = df['binding_affinity'] * -1
    df['lig_rmsd'] = df['lig_rmsd'] * -1
    # Standardize the data
    df.iloc[:, 3:] = (df.iloc[:, 3:] - df.iloc[:, 3:].mean()) / df.iloc[:, 3:].std()
    
    # Create a correlation plot to check for multicollinearity
    print('Generating correlation plot...')
    print('\nClose the plot window to continue. The plot will be saved as "correlation_plot.png".')
    corr = df.iloc[:, 3:].corr()
    sns.set_theme(style="white")
    sns.set_style("whitegrid")
    sns.heatmap(corr, xticklabels=corr.columns, yticklabels=corr.columns, cmap='Greens')
    # Display the values of the correlation coefficients
    for i in range(len(corr)):
        for j in range(len(corr)):
            text = plt.text(j + 0.5, i + 0.5, round(corr.iloc[i, j], 2), ha="center", va="center", color="black")
    plt.title('Correlation Plot')
    plt.tight_layout()
    # Save the figure
    plt.savefig(file_dir + '\correlation_plot.png', dpi=300)
    plt.show()
    
    # Specify the k-means parameters
    kmeans_kwargs = {
        "init": "k-means++",
        "n_init": "auto",
        "random_state": 0,
    }
    
    # Generate an elbow plot to determine the number of clusters
    print('Generating elbow plot...')
    print('\nClose the plot window to continue. The plot will be saved as "elbow_plot.png".')
    sse = []
    for k in range(1, 11):
        kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
        kmeans.fit(df.iloc[:, 3:])
        sse.append(kmeans.inertia_)
    plt.style.use("fivethirtyeight")
    plt.rcParams['axes.facecolor'] = 'white'
    plt.plot(range(1, 11), sse, 'bx-')
    plt.title('Elbow Plot')
    plt.xticks(range(1, 11))
    plt.xlabel("Number of Clusters")
    plt.ylabel("SSE")
    plt.tight_layout()
    # Save the figure
    plt.savefig(file_dir + '\elbow_plot.png', dpi=300)
    plt.show()
    # Determine the optimal number of clusters programmatically
    kl = KneeLocator(range(1, 11), sse, curve="convex", direction="decreasing")
    # Print the optimal number of clusters
    print('\nAccording to the elbow plot, the optimal number of clusters is: ' + str(kl.elbow))
    
    # Generate a silhouette plot to determine the number of clusters
    print('Generating silhouette plot...')
    print('\nClose the plot window to continue. The plot will be saved as "silhouette_plot.png".')
    silhouette_scores = []
    for k in range(2, 11):
        kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
        kmeans.fit(df.iloc[:, 3:])
        score = silhouette_score(df.iloc[:, 3:], kmeans.labels_)
        silhouette_scores.append(score)
    # Plot the silhouette scores
    plt.style.use("fivethirtyeight")
    plt.rcParams['axes.facecolor'] = 'white'
    plt.plot(range(2, 11), silhouette_scores, 'bx-')
    plt.title('Silhouette Plot')
    plt.xticks(range(2, 11))
    plt.xlabel("Number of Clusters")
    plt.ylabel("Silhouette Score")
    plt.tight_layout()
    # Save the figure
    plt.savefig(file_dir + '\silhouette_plot.png', dpi=300)
    plt.show()
    # Determine the cluster with the highest silhouette score
    max_score = max(silhouette_scores)
    max_score_index = silhouette_scores.index(max_score)
    # Print the optimal number of clusters
    print('\nAccording to the silhouette plot, the optimal number of clusters is: ' + str(max_score_index + 2))
    
    # Specify the k-means parameters for the final clustering
    n_clusters = input("Enter the number of clusters: ")
    n_clusters = int(n_clusters)
    kmeans = KMeans(n_clusters=n_clusters, **kmeans_kwargs)
    kmeans.fit(df.iloc[:, 3:])
    # Get the cluster labels
    labels = kmeans.labels_
    
    # Perform PCA on the standardized data
    print('\nPerforming PCA to reduce the data to 2 dimensions...')
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(df.iloc[:, 3:])
    # Determine the loadings of the principal components
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    # Print the loadings of the principal components
    print('\nThe loadings of the principal components are:\n')
    print(pd.DataFrame(loadings, columns=['PC1', 'PC2'], index=df.iloc[:, 3:].columns))
    pd.DataFrame(loadings, columns=['PC1', 'PC2'], index=df.iloc[:, 3:].columns).to_csv(file_dir + '\pca_loadings.csv')
    # Add the principal components to the dataframe
    df['PCA1'] = principal_components[:, 0]
    df['PCA2'] = principal_components[:, 1]
    # Add the cluster labels to the dataframe
    df['cluster'] = labels
    # Write the new dataframe to a csv file
    print('\nWriting the new dataframe to a csv file.\nThe file will be saved as "kmeans_clustering_summary.csv".')
    df.to_csv(file_dir + '\kmeans_clustering_summary.csv', index=False)
    
    # Ask the user for the name of the self-docking pose
    self_dock_species = input("\nEnter the species name of the self-docking pose to be used as a pseudo-reference (case sensitive): ")
    while self_dock_species not in df['species'].unique().tolist():
        self_dock_species = input("The species name entered is not found in the dataframe. Please enter the species name again: ")
    # Ask the user for the pose number of the self-docking pose
    self_dock_pose = input("Enter the pose number of the self-docking pose to be used as a pseudo-reference: ")
    df.columns = [col.upper() for col in df.columns]

    # Plot the clusters in 2D using seaborn
    print('\nPlotting the clusters in 2D...')
    # Create a cmap for the clusters
    cmap = ListedColormap(['#ff7f0e', '#2ca02c', "#1f77b4"])
    sns.set(font_scale=1.2)
    sns.set_style('whitegrid')
    sns.scatterplot(x='PCA1', y='PCA2', hue='CLUSTER', data=df, palette=cmap, legend='full', s=100)
    plt.xlabel('PC1 ({:.2f}%)'.format(pca.explained_variance_ratio_[0]*100))
    plt.ylabel('PC2 ({:.2f}%)'.format(pca.explained_variance_ratio_[1]*100))
    plt.scatter(df.loc[(df['SPECIES'] == self_dock_species) & (df['POSE'] == int(self_dock_pose)), 'PCA1'], df.loc[(df['SPECIES'] == self_dock_species) & (df['POSE'] == int(self_dock_pose)), 'PCA2'], s=200, facecolors='none', edgecolors='black', linewidths=2)
    plt.tight_layout()
    # Save the figure
    plt.savefig(file_dir + '\\training_set_clustering.png', dpi=300)
    plt.show()
    
    # Find the cluster for the self-docking pose
    ref_cluster = df.loc[(df['SPECIES'] == self_dock_species) & (df['POSE'] == int(self_dock_pose)), 'CLUSTER'].values[0]
    # Filter the dataframe to only include the reference cluster
    df0 = df[df['CLUSTER'] == ref_cluster]
    # Get the list of species in cluster 0
    species0 = df0['SPECIES'].unique().tolist()
    # Print as susceptible species
    print('Susceptible species:')
    for i in species0:
        print(i)
    # Determine the species that are not in cluster 0
    species1 = df['SPECIES'].unique().tolist()
    species1 = [x for x in species1 if x not in species0]
    # Print as non susceptible species
    print('Non-susceptible species:')
    for i in species1:
        print(i)
    # Write the susceptible species summary to a text file
    with open(file_dir + '\susceptibility_summary.txt', 'w') as f:
        f.write('Susceptible species:\n')
        for i in species0:
            f.write(i + '\n')
        f.write('\nNon-susceptible species:\n')
        for i in species1:
            f.write(i + '\n')
    
    print("\nAll of the plots for this analysis have been saved to {}".format(file_dir))
    print("A summary of the species susceptibility calls has been saved as 'susceptibility_summary.txt'.")



    # Ask the user if they would like to make predictions on a test set
    test_set = input("\nWould you like to make predictions on a test set? (y/n): ")
    test_set = test_set.lower()
    while test_set not in ['y', 'n']:
        test_set = input("Please enter 'y' or 'n': ")
    if test_set == 'n':
        return
    # Specify the features as the 4 docking metrics, prior to PCA
    X = df.iloc[:, 3:7]
    # Use the cluster labels as the target
    y = df['CLUSTER']
    while True:
        # Split the data into training and testing sets to train the SVM model
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=109)
        X_val, X_test, y_val, y_test = train_test_split(X_test, y_test, test_size=0.5, random_state=109)
        # Print the shapes of the training, validation, and testing sets
        print('\nTraining set shape: {}'.format(X_train.shape))
        print('Validation set shape: {}'.format(X_val.shape))
        print('Testing set shape: {}'.format(X_test.shape))
        # Determine the optimal kernel for the SVM model
        print('\nDetermining the optimal kernel for the SVM model...')
        kernels = ['linear', 'poly', 'rbf', 'sigmoid']
        for kernel in kernels:
            clf = svm.SVC(kernel=kernel)
            clf.fit(X_train, y_train)
            y_pred = clf.predict(X_test)
            print('\nKernel: {}'.format(kernel))
            print("Model accuracy:", metrics.accuracy_score(y_test, y_pred))
            print("Model precision:", metrics.precision_score(y_test, y_pred, average='weighted'))
            print("Model recall:", metrics.recall_score(y_test, y_pred, average='weighted'))
            print("Model F1 score:", metrics.f1_score(y_test, y_pred, average='weighted'))
        # Ask the user for the optimal kernel
        kern = input("\nEnter the optimal kernel for the SVM model: ")
        while kern not in kernels:
            kernel = input("Please enter a valid kernel: ")
        # Evaluate the validation set using the optimal kernel
        print('\nEvaluating the validation set using the optimal kernel...')
        clf = svm.SVC(kernel=kern)
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_val)
        print('\nKernel: {}'.format(kern))
        print("Model accuracy:", metrics.accuracy_score(y_val, y_pred))
        print("Model precision:", metrics.precision_score(y_val, y_pred, average='weighted'))
        print("Model recall:", metrics.recall_score(y_val, y_pred, average='weighted'))
        print("Model F1 score:", metrics.f1_score(y_val, y_pred, average='weighted'))
        model_sele = input("\n Based on this evaluation, would you like to use this model to make predictions on a test set? (y/n): ")
        model_sele = model_sele.lower()
        while model_sele not in ['y', 'n']:
            model_sele = input("Please enter 'y' or 'n': ")
        if model_sele == 'y':
            break
        elif model_sele == 'n':
            print('\nPredictions will not be made on a test set.')
            return
    # Create a SVM Classifier
    clf = svm.SVC(kernel=kern)
    clf = clf.fit(X_train, y_train)
    # Ask for the file path to a test set to make predictions on
    test_set = input("\nEnter the file path to a test set to make predictions on: ")
    test_set = test_set.replace('"', "")
    test_df = pd.read_csv(test_set)
    # Prep and standardize the test set
    test_df['binding_affinity'] = test_df['binding_affinity'] * -1
    test_df['lig_rmsd'] = test_df['lig_rmsd'] * -1
    test_df.iloc[:, 3:] = (test_df.iloc[:, 3:] - test_df.iloc[:, 3:].mean()) / test_df.iloc[:, 3:].std()
    # Make predictions on the test set using the SVM model and the PCA data
    test_df.columns = [col.upper() for col in test_df.columns]
    y_pred = clf.predict(test_df.iloc[:, 3:7])
    # Add the predictions to the dataframe
    test_df['CLUSTER'] = y_pred

    # Remove the PCA columns by reloading the summary file
    df = pd.read_csv(summary_file)
    # Multiply the "binding_affinity" and "lig_rmsd" columns by -1 so that higher values indicate a better pose
    df['binding_affinity'] = df['binding_affinity'] * -1
    df['lig_rmsd'] = df['lig_rmsd'] * -1
    # Standardize the data
    df.iloc[:, 3:] = (df.iloc[:, 3:] - df.iloc[:, 3:].mean()) / df.iloc[:, 3:].std()
    # Add the cluster labels to the dataframe
    df['CLUSTER'] = y
    # Convert all column names to uppercase
    df.columns = [col.upper() for col in df.columns]
    
    # Combine with the original dataframe
    df = pd.concat([df, test_df])
    # Perform PCA on the combined dataframe for visualization
    print('\nPerforming PCA to reduce the data to 3 dimensions...')
    pca = PCA(n_components=3)
    principal_components = pca.fit_transform(df.iloc[:, 3:7])
    # Add the principal components to the dataframe
    df['PCA1'] = principal_components[:, 0]
    df['PCA2'] = principal_components[:, 1]
    df['PCA3'] = principal_components[:, 2]

    # Determine the loadings of the principal components
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    # Print the loadings of the principal components
    print('\nThe loadings of the principal components are:\n')
    print(pd.DataFrame(loadings, columns=['PC1', 'PC2', 'PC3'], index=df.iloc[:, 3:7].columns))
    pd.DataFrame(loadings, columns=['PC1 ({:.2f}%)'.format(pca.explained_variance_ratio_[0]*100), 'PC2 ({:.2f}%)'.format(pca.explained_variance_ratio_[1]*100), 'PC3 ({:.2f}%)'.format(pca.explained_variance_ratio_[2]*100)], index=df.iloc[:, 3:7].columns).to_csv(file_dir + '\pca_loadings.csv')

    # Plot the combined data in 2D
    print('\nPlotting the clusters in 2D...')
    print('\nClose the plot window to continue. The plot will be saved as "combined_cluster_plot.png".')
    sns.set(font_scale=1.2)
    sns.set_style("whitegrid")
    sns.scatterplot(x='PCA1', y='PCA2', hue='CLUSTER', data=df, palette=cmap, s=100)
    plt.xlabel('PC1 ({:.2f}%)'.format(pca.explained_variance_ratio_[0]*100))
    plt.ylabel('PC2 ({:.2f}%)'.format(pca.explained_variance_ratio_[1]*100))
    # Label the added test set
    sns.scatterplot(x='PCA1', y='PCA2', data=df.iloc[-len(test_df):, :], s=200, marker='o', facecolors='none', edgecolor='black', linewidth=2)
    plt.tight_layout()
    # Save the figure
    plt.savefig(file_dir + '\\predicted_clustering_2D.png', dpi=300)
    plt.show()

    # Generate a 3D scatter plot of the combined data
    plot_3d = input("\nWould you like to plot the clusters in 3D? (y/n): ")
    while plot_3d not in ['y', 'n']:
        plot_3d = input("Please enter a valid response (y/n): ")
    if plot_3d == 'y':
        print('\nPlotting the clusters in 3D...')
        print('\nClose the plot window to continue. The plot will be saved as "combined_cluster_plot_3D.png".')
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(df['PCA1'], df['PCA2'], df['PCA3'], c=df['CLUSTER'], cmap=cmap, s=100, alpha=0.6)
        ax.set_xlabel('PC1 ({:.2f}%)'.format(pca.explained_variance_ratio_[0]*100))
        ax.set_ylabel('PC2 ({:.2f}%)'.format(pca.explained_variance_ratio_[1]*100))
        ax.set_zlabel('PC3 ({:.2f}%)'.format(pca.explained_variance_ratio_[2]*100))
        # Label the added test set
        ax.scatter(df.iloc[-len(test_df):, -3], df.iloc[-len(test_df):, -2], df.iloc[-len(test_df):, -1], s=200, marker='o', facecolors='none', edgecolor='black', linewidth=2)
        plt.tight_layout()
        # Save the figure
        plt.savefig(file_dir + '\\predicted_clustering_3D', dpi=300)
        plt.show()

    # Find the cluster for the self-docking pose
    ref_cluster = df.loc[(df['SPECIES'] == self_dock_species) & (df['POSE'] == int(self_dock_pose)), 'CLUSTER'].values[0]
    # Filter the dataframe to only include the reference cluster
    df0 = df[df['CLUSTER'] == ref_cluster]
    # Get the list of species in cluster 0
    species0 = df0['SPECIES'].unique().tolist()
    # Only include species from the test set
    species1 = [x for x in species0 if x in test_df['SPECIES'].unique().tolist()]
    print('Susceptible species:')
    for i in species1:
        print(i)
    # Determine the species that are not in cluster 0
    species2 = df['SPECIES'].unique().tolist()
    species2 = [x for x in species1 if x not in species0]
    # Print as non susceptible species
    print('Non-susceptible species:')
    for i in species2:
        print(i)
    # Write the susceptible species summary to a text file
    with open(file_dir + '\susceptibility_summary.txt', 'w') as f:
        f.write('Susceptible species:\n')
        for i in species1:
            f.write(i + '\n')
        f.write('\nNon-susceptible species:\n')
        for i in species2:
            f.write(i + '\n')
    print("\nAll of the plots for this analysis have been saved to {}".format(file_dir))
    print("A summary of the species susceptibility calls has been saved as 'susceptibility_summary.txt'.")

cluster_analysis()