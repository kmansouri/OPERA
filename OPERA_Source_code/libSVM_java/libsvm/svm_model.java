package libsvm;

// Class representing the SVM model
public class svm_model implements java.io.Serializable {
    
    public svm_parameter param;  
    // SVM parameters used to train the model
    
    // For classification
    public int nr_class;        
    // The number of classes in the dataset.
    // For binary classification, nr_class = 2. For multiclass classification, it is greater than 2.
    
    public int l;               
    // The total number of Support Vectors (SV) in the model.
    
    public svm_node[][] SV;     
    // Support Vectors (SV) of the model.
    // SV is a 2-dimensional array where each element represents a data point (svm_node) in the SV set.
    
    public double[][] sv_coef;  
    // Coefficients for Support Vectors in decision functions.
    // sv_coef is a 2-dimensional array where sv_coef[k-1][l] represents the coefficient of the l-th Support Vector for the k-th class.
    // For binary classification, there's only one coefficient (k=1) for each Support Vector.
    
    public double[] rho;        
    // Constants in decision functions.
    // rho is an array containing the constants used in the decision functions.
    // For binary classification, there's only one constant.
    // For multiclass classification, there are (nr_class * (nr_class-1)) / 2 constants.
    // These constants are used to calculate the decision values for classification.
    
    public double[] probA;         
    // Pairwise probability information used in probability estimates.
    // probA is an array containing information used to estimate class membership probabilities in classification.
    
    public double[] probB;
    // Pairwise probability information used in probability estimates.
    // probB is another array containing information used to estimate class membership probabilities in classification.
    
    public int[] sv_indices;       
    // Indices of Support Vectors (SV) in the training set.
    // sv_indices is an array containing values in the range [1, ..., num_training_data] to indicate the positions of Support Vectors in the original training data.
    
    // For classification only
    
    public int[] label;          
    // Labels of each class.
    // label is an array where label[k] represents the label of the k-th class.
    
    public int[] nSV;           
    // Number of Support Vectors (SV) for each class.
    // nSV is an array where nSV[k] represents the number of Support Vectors for the k-th class.
    // The sum of all elements in nSV should be equal to 'l', the total number of SVs.
};
