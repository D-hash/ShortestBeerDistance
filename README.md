# ShortestBeerDistance
This is the implementation of BUILD-WHL algorithm presented in "Scalable Graph Indexing for Faster Shortest Beer Path Queries" submitted at Very Large DataBase Conference (VLDB) 2025.
Given a weighted beer graph G=(V,E), it computes a weighted Highway Labeling of G that is exploited at runtime to answer beer distance queries by an ad-hoc query algorithm.

The code extends the implementation of Farhan https://github.com/mufarhan/highway_labelling

The current implementation requires Networkit library https://github.com/networkit/networkit
