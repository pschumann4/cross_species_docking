"""
 This script will perform a k-nearest neighbors cluster analysis 
 for generating species susceptibility calls based on the docking 
 metrics of the self-docking pose (i.e., reference pose). The user 
 will need to provide a summary file containing the docking      
 metrics for the self-docking pose and the poses of the test set. 
 The user will also need to provide the species name and pose number 
 of the self-docking pose. This script will output a summary of the 
 species susceptibility calls and plots of the clustering results.                                                                              
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from kneed import KneeLocator
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.neighbors import KNeighborsClassifier


def cluster_analysis():
    """
    This is the main function that will be called by the user.
    It will prompt the user for the file directory containing the summary file,
    the species name and pose number of the self-docking pose, the number of clusters,
    and if they want to perform a k-nearest neighbors analysis.
    """
    # Read the summary csv file
    summary_file = input("Enter the file path for the ensemble-docking summary file: ")
    # Remove the quotes from the file path
    summary_file = summary_file.replace('"', "")
    file_dir = os.path.dirname(summary_file)
    df = pd.read_csv(summary_file)
    # Multiply the "binding_affinity" and "lig_rmsd" columns by -1
    # so that higher values indicate a better pose
    df["binding_affinity"] = df["binding_affinity"] * -1
    df["lig_rmsd"] = df["lig_rmsd"] * -1
    # Standardize the data
    df.iloc[:, 3:] = (df.iloc[:, 3:] - df.iloc[:, 3:].mean()) / df.iloc[:, 3:].std()

    # Create a correlation plot to check for multicollinearity
    print("Generating correlation plot...")
    print(
        "\nClose the plot window to continue. "
        'The plot will be saved as "correlation_plot.png".'
    )
    corr = df.iloc[:, 3:].corr()
    sns.set_theme(style="white")
    sns.set_style("whitegrid")
    sns.heatmap(corr, xticklabels=corr.columns, yticklabels=corr.columns, cmap="Greens")
    # Display the values of the correlation coefficients
    for i in range(len(corr)):
        for j in range(len(corr)):
            plt.text(
                j + 0.5,
                i + 0.5,
                round(corr.iloc[i, j], 2),
                ha="center",
                va="center",
                color="black",
            )
    # Change the label plif_tanimoto to PLIF Tanimoto
    labels = [item.get_text() for item in plt.gca().get_xticklabels()]
    labels = ["Binding Affinity", "PPS-Score", "Ligand RMSD", "PLIF Tc"]
    plt.gca().set_xticklabels(labels)
    plt.gca().set_yticklabels(labels)
    plt.tight_layout()
    # Save the figure
    plt.savefig(file_dir + "\correlation_plot.png", dpi=300)
    plt.show()

    # Specify the k-means parameters
    kmeans_kwargs = {
        "init": "k-means++",
        "n_init": "auto",
        "random_state": 40,
    }

    # Generate an elbow plot to determine the number of clusters
    print("Generating elbow plot...")
    print(
        "\nClose the plot window to continue. "
        'The plot will be saved as "elbow_plot.png".'
    )
    sse = []
    for k in range(1, 11):
        kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
        kmeans.fit(df.iloc[:, 3:])
        sse.append(kmeans.inertia_)
    # Determine the optimal number of clusters programmatically
    kl = KneeLocator(range(1, 11), sse, curve="convex", direction="decreasing")
    plt.style.use("seaborn-v0_8-colorblind")
    plt.plot(range(1, 11), sse, "bx-")
    plt.xticks(range(1, 11))
    plt.xlabel("Number of Clusters")
    plt.ylabel("SSE")
    plt.gca().lines[0].set_color("black")
    plt.gca().lines[0].set_linewidth(2)
    # Plot a red dot at the optimal number of clusters
    plt.plot(
        kl.elbow,
        sse[kl.elbow - 1],
        marker="o",
        markersize=10,
        markerfacecolor="red",
        markeredgecolor="black"
    )
    plt.tight_layout()
    # Save the figure
    plt.savefig(file_dir + "\elbow_plot.png", dpi=300)
    plt.show()
    
    # Print the optimal number of clusters
    print(
        "\nAccording to the elbow plot, "
        "the optimal number of clusters is: " + str(kl.elbow)
    )

    # Generate a silhouette plot to determine the number of clusters
    print("Generating silhouette plot...")
    print(
        "\nClose the plot window to continue. "
        'The plot will be saved as "silhouette_plot.png".'
    )
    silhouette_scores = []
    for k in range(2, 11):
        kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
        kmeans.fit(df.iloc[:, 3:])
        score = silhouette_score(df.iloc[:, 3:], kmeans.labels_)
        silhouette_scores.append(score)
    # Determine the cluster with the highest silhouette score
    max_score = max(silhouette_scores)
    max_score_index = silhouette_scores.index(max_score)
    # Plot the silhouette scores
    plt.style.use("seaborn-v0_8-colorblind")
    plt.plot(range(2, 11), silhouette_scores, "bx-")
    plt.xticks(range(2, 11))
    plt.xlabel("Number of Clusters")
    plt.ylabel("Silhouette Score")
    plt.gca().lines[0].set_color("black")
    plt.gca().lines[0].set_linewidth(2)
    # Plot a red dot at the optimal number of clusters
    plt.plot(
        max_score_index + 2,
        silhouette_scores[max_score_index],
        marker="o",
        markersize=10,
        markerfacecolor="red",
        markeredgecolor="black"
    )
    plt.tight_layout()
    # Save the figure
    plt.savefig(file_dir + os.sep + "silhouette_plot.png", dpi=300)
    plt.show()

    # Print the optimal number of clusters
    print(
        "\nAccording to the silhouette plot, "
        "the optimal number of clusters is: " + str(max_score_index + 2)
    )

    # Specify the k-means parameters for the final clustering
    n_clusters = input("Enter the number of clusters: ")
    n_clusters = int(n_clusters)
    kmeans = KMeans(n_clusters=n_clusters, **kmeans_kwargs)
    kmeans.fit(df.iloc[:, 3:])
    # Get the cluster labels
    labels = kmeans.labels_

    # Perform PCA on the standardized data
    print("\nPerforming PCA to reduce the data to 2 dimensions...")
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(df.iloc[:, 3:])

    # Add the principal components to the dataframe
    df["PCA1"] = principal_components[:, 0]
    df["PCA2"] = principal_components[:, 1]
    # Add the cluster labels to the dataframe
    df["cluster"] = labels
    # Write the new dataframe to a csv file
    print(
        "\nWriting the new dataframe to a csv file."
        '\nThe file will be saved as "kmeans_clustering_summary.csv".'
    )
    df.to_csv(file_dir + "\kmeans_clustering_summary.csv", index=False)

    # Ask the user for the name of the self-docking pose
    self_dock_species = input("Enter the species name of the self-docking pose to be used as a pseudo-reference (case sensitive): ")
    while self_dock_species not in df["species"].unique().tolist():
        self_dock_species = input("The species name entered is not found in the dataframe. Please enter the species name again: ")
    # Ask the user for the pose number of the self-docking pose
    self_dock_pose = input("Enter the pose number of the self-docking pose to be used as a pseudo-reference: ")
    df.columns = [col.upper() for col in df.columns]

    # Plot the clusters in 2D using seaborn
    print("\nPlotting the clusters in 2D...")
    # Create a cmap for the clusters
    cmap = ListedColormap(["#D81B60", "#1E88E5", "#FFC107", "#004D40", "#CC79A7", "#56B4E9"])
    sns.set(font_scale=1.2)
    sns.set_style("whitegrid")
    sns.scatterplot(
        x="PCA1", y="PCA2", hue="CLUSTER", data=df, palette=cmap.colors, legend="full", s=75
    )
    sns.move_legend(plt.gca(), "upper center", bbox_to_anchor=(0.5, 1.15), ncol = len(df["CLUSTER"].unique().tolist()), title=None)
    plt.xlabel("PC1 ({:.2f}%)".format(pca.explained_variance_ratio_[0] * 100))
    plt.ylabel("PC2 ({:.2f}%)".format(pca.explained_variance_ratio_[1] * 100))
    plt.scatter(
        df.loc[
            (df["SPECIES"] == self_dock_species) & (df["POSE"] == int(self_dock_pose)),
            "PCA1",
        ],
        df.loc[
            (df["SPECIES"] == self_dock_species) & (df["POSE"] == int(self_dock_pose)),
            "PCA2",
        ],
        s=100,
        facecolors="none",
        edgecolors="black",
        linewidths=2,
    )
    plt.tight_layout()
    # Save the figure
    plt.savefig(file_dir + "\\training_set_clustering.png", dpi=300)
    plt.show()

    print("\nAll of the plots for this analysis have been saved to {}".format(file_dir))

    # Ask the user if they would like to make predictions on a test set
    test = input("\nWould you like to make predictions on a test set? (y/n): ")
    test = test.lower()
    while test not in ["y", "n"]:
        test = input("Please enter 'y' or 'n': ")
    if test == "n":
        return

    # Ask for the file path to a test set to make predictions on
    test_set = input("\nEnter the file path to a test set to make predictions on: ")
    test_set = test_set.replace('"', "")
    file_dir = os.path.dirname(test_set)
    test_df = pd.read_csv(test_set)
    # Prep and standardize the test set
    test_df["binding_affinity"] = test_df["binding_affinity"] * -1
    test_df["lig_rmsd"] = test_df["lig_rmsd"] * -1
    test_df.iloc[:, 3:] = (
        test_df.iloc[:, 3:] - test_df.iloc[:, 3:].mean()
    ) / test_df.iloc[:, 3:].std()
    # Convert all column names to uppercase
    test_df.columns = [col.upper() for col in test_df.columns]

    # Specify the features as the 4 docking metrics, prior to PCA
    X = df.iloc[:, 3:7]
    # Use the cluster labels as the target
    Y = df["CLUSTER"]
    # Use k-nearest neighbors to predict the cluster labels for the test set
    print("\nDetermining optimal number of neighbors via 5-fold cross-validation...")

    # Determine the optimal number of neighbors via 5-fold cross-validation
    k_range = list(range(1, round(len(df) / 2)))
    # Create a dictionary of the parameter values to test
    param_grid = dict(n_neighbors=k_range)
    # Instantiate the grid
    grid = GridSearchCV(
        KNeighborsClassifier(),
        param_grid,
        cv=KFold(n_splits=5, shuffle=True, random_state=42),
        scoring="accuracy"
    )
    # Fit the grid with data
    grid.fit(X, Y)
    # Plot the results
    plt.style.use("seaborn-v0_8-colorblind")
    plt.plot(k_range, grid.cv_results_["mean_test_score"])
    plt.gca().lines[0].set_color("black")
    plt.gca().lines[0].set_linewidth(2)
    plt.xlabel("Value of K")
    plt.ylabel("Cross-Validated Accuracy")
    plt.tight_layout()
    plt.savefig(file_dir + "\\knn_cross_validation.png", dpi=300)
    plt.show()

    # View the results
    print(
        "The optimal number of neighbors is: {}".format(
            grid.best_params_["n_neighbors"]
        )
    )

    # Write the cv.results_ to a csv file
    cv_results = pd.DataFrame(grid.cv_results_)
    cv_results.to_csv(file_dir + "\\knn_validation_results.csv")

    # Instantiate the model with the optimal number of neighbors
    knn = KNeighborsClassifier(n_neighbors=grid.best_params_["n_neighbors"])
    # Fit the model
    knn.fit(X, Y)

    # Make predictions on the test set
    print("\nMaking predictions on the test set...")
    Z = knn.predict(test_df.iloc[:, 3:7])
    # Add the cluster labels to the dataframe
    test_df["CLUSTER"] = Z

    # Remove the PCA columns by reloading the summary file
    train_df = pd.read_csv(summary_file)
    train_df["binding_affinity"] = train_df["binding_affinity"] * -1
    train_df["lig_rmsd"] = train_df["lig_rmsd"] * -1
    train_df.iloc[:, 3:] = (train_df.iloc[:, 3:] - train_df.iloc[:, 3:].mean()) / train_df.iloc[:, 3:].std()
    train_df.columns = [col.upper() for col in train_df.columns]
    # Add the cluster labels
    train_df["CLUSTER"] = Y

    # Combine with the original dataframe
    df = pd.concat([train_df, test_df])

    # Perform PCA on the combined dataframe for visualization
    print("\nPerforming PCA to reduce the data to 3 dimensions...")
    pca = PCA(n_components=3)
    principal_components = pca.fit_transform(df.iloc[:, 3:7])
    # Add the principal components to the dataframe
    df["PCA1"] = principal_components[:, 0]
    df["PCA2"] = principal_components[:, 1]
    df["PCA3"] = principal_components[:, 2]

    # Determine the loadings of the principal components
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    # Print the loadings of the principal components
    print("\nThe loadings of the principal components are:\n")
    print(
        pd.DataFrame(
            loadings, columns=["PC1", "PC2", "PC3"], index=df.iloc[:, 3:7].columns
        )
    )
    pd.DataFrame(
        loadings,
        columns=[
            "PC1 ({:.2f}%)".format(pca.explained_variance_ratio_[0] * 100),
            "PC2 ({:.2f}%)".format(pca.explained_variance_ratio_[1] * 100),
            "PC3 ({:.2f}%)".format(pca.explained_variance_ratio_[2] * 100),
        ],
        index=df.iloc[:, 3:7].columns,
    ).to_csv(file_dir + "\pca_loadings.csv")

    # Save the dataframe to a csv file
    df.to_csv(file_dir + "\\predicted_clustering_summary.csv", index=False)

    # Get a list of all the species names from the train_df
    train_species = train_df["SPECIES"].unique().tolist()

    # Plot the combined data in 2D
    print("\nPlotting the clusters in 2D...")
    print(
        '\nClose the plot window to continue. The plot will be saved as "combined_cluster_plot.png".'
    )
    sns.set(font_scale=1.2)
    sns.set_style("whitegrid")
    sns.scatterplot(
        x="PCA1", y="PCA2", hue="CLUSTER", data=df, palette=cmap.colors, legend="full", s=35
    )
    sns.move_legend(plt.gca(), "upper center", bbox_to_anchor=(0.5, 1.15), ncol = len(df["CLUSTER"].unique().tolist()), title=None)
    plt.xlabel("PC1 ({:.2f}%)".format(pca.explained_variance_ratio_[0] * 100))
    plt.ylabel("PC2 ({:.2f}%)".format(pca.explained_variance_ratio_[1] * 100))
    # Label all the species in found in the train_species list
    for i in train_species:
        plt.scatter(
            df.loc[(df["SPECIES"] == i), "PCA1"],
            df.loc[(df["SPECIES"] == i), "PCA2"],
            s=35,
            facecolors="none",
            edgecolors="black",
            linewidths=1,
        )
    plt.tight_layout()
    # Save the figure
    plt.savefig(file_dir + '\\predicted_clustering_2D.png', dpi=300)
    plt.show()

    # Generate a 3D scatter plot of the combined data
    plot_3d = input("\nWould you like to plot the clusters in 3D? (y/n): ")
    while plot_3d not in ["y", "n"]:
        plot_3d = input("Please enter a valid response (y/n): ")
    if plot_3d == "y":
        print("\nPlotting the clusters in 3D...")
        print(
            "\nClose the plot window to continue. "
            'The plot will be saved as "combined_cluster_plot_3D.png".'
        )
        cmap = LinearSegmentedColormap.from_list(
            "custom", cmap.colors[: len(df["CLUSTER"].unique().tolist())]
        )
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(
            df["PCA1"],
            df["PCA2"],
            df["PCA3"],
            c=df["CLUSTER"],
            cmap=cmap,
            s=35,
            alpha=0.65,
        )
        # Label all the species in found in the train_species list
        for i in train_species:
            ax.scatter(
                df.loc[(df["SPECIES"] == i), "PCA1"],
                df.loc[(df["SPECIES"] == i), "PCA2"],
                df.loc[(df["SPECIES"] == i), "PCA3"],
                s=35,
                alpha=1,
                facecolors="none",
                edgecolors="black",
                linewidths=2,
            )
        ax.set_xlabel("PC1 ({:.2f}%)".format(pca.explained_variance_ratio_[0] * 100))
        ax.set_ylabel("PC2 ({:.2f}%)".format(pca.explained_variance_ratio_[1] * 100))
        ax.set_zlabel("PC3 ({:.2f}%)".format(pca.explained_variance_ratio_[2] * 100))
        plt.tight_layout()
        # Save the figure
        plt.savefig(file_dir + "\\predicted_clustering_3D", dpi=300)
        plt.show()

    # Find the cluster for the self-docking pose
    ref_cluster = df.loc[
        (df["SPECIES"] == self_dock_species) & (df["POSE"] == int(self_dock_pose)),
        "CLUSTER",
    ].values[0]
    # Filter the dataframe to only include the reference cluster
    df_ref = df[df["CLUSTER"] == ref_cluster]
    # Get the list of species in cluster 0
    sus_species = df_ref["SPECIES"].unique().tolist()
    print("Susceptible species:")
    for i in sus_species:
        print(i)
    # Determine the species that are not in cluster 0
    nsus_species = df["SPECIES"].unique().tolist()
    # Remove the species in species_ref from species_nsus
    nsus_species = [x for x in nsus_species if x not in sus_species]
    # Print as non susceptible species
    print("Non-susceptible species:")
    for i in nsus_species:
        print(i)
    # Write the susceptible species summary to a text file
    with open(file_dir + "\susceptibility_summary.txt", "w") as f:
        f.write("Susceptible species:\n")
        for i in sus_species:
            f.write(i + "\n")
        f.write("\nNon-susceptible species:\n")
        for i in nsus_species:
            f.write(i + "\n")
    print(
        "\nAll of the plots for this analysis have " "been saved to {}".format(file_dir)
    )
    print(
        "A summary of the species susceptibility calls has been "
        "saved as 'susceptibility_summary.txt'."
    )

if __name__ == "__main__":
    cluster_analysis()
