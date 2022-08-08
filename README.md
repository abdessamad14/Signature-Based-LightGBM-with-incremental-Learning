Signature-Based LightGBM with incremental Learning
==============================

Today, real-time data streams receive more and more attention in many areas. Storing and then learning data in this context can lead to several technical difficulties. Firstly, the storage space will be larger and larger over time. Second, it is necessary to train the entire machine learning model from the beginning as new batches of data come in and regularly reconstruct new models from scratch. Incremental Learning allows us to use information as soon as it is available leading to all-time up-to-date models and also reducing the costs for data storage and maintenance. 

In this master's thesis, a Signature-Based LightGBM with incremental learning is proposed. This is a method that detects at the first place concept drift on Time Series Data Stream in a robust way. Our drift detector monitors the drift incrementally based on confidence scores of LightGBM model and not the data which might not be available during the production phase. Once a drift is detected, the new knowledge is incrementally learned using the LightGBM Classifier. Using Signature Transformation of Time Series Data Stream as well as Incremental Bayesian Hyper-parameter Optimization is our model resilience against Catastrophic forgetting. 

