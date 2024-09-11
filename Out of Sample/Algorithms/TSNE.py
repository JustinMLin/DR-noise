from openTSNE import TSNE

def tsne(Z, perplexity=30, n_components=2, exaggeration=None):
    return TSNE(perplexity=perplexity,
                n_components=n_components,
                initialization="pca",
                exaggeration=exaggeration,
                negative_gradient_method="auto").fit(Z)
    
def predict_tsne(tsne, new_data):
    return tsne.transform(new_data)