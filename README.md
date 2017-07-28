# sharing methods and interesting science

CCI.py is a command-line demo of a tool I use to "clean up" outlier streamlines. Each streamline gets a vote in each other streamline's "cluster confidence index" score based on the similarity of the two paths. Streamlines that follow paths with many other streamlines have higher confidence than those following a pathway alone, in the context of the dataset of streamlines. The provided code packages the CCI score into the .trk file so that the streamline dataset can be opened in trackvis and filtered by CCI with a slider bar.

pathlength_batch.py is a demo of how to run the path_length function added to Diffusion Imaging in Python to generate path-length niftis from a set of streamlines saved in a .trk file and a nifti with the binary ROI from which the pathlength should be measured. This code is being used for a Radiation Therapy project at UCSF trying to predict the recurrence of brain tumors based on white matter pathways, modeled by tractography (https://github.com/nipy/dipy/pull/1114).

neurosurgery_trees is a fun demo re-analyzing a study (http://thejns.org/doi/abs/10.3171/2015.6.JNS142203) predicting post-surgical deficits based on pre- and post-surgical tractography modeling of relevant white matter structures using decision trees. The notebook (Deficit_Prediction_EduCohort.ipyb) summarizes results.
