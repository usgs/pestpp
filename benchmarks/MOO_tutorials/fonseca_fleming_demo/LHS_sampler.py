import numpy as np
import os
from scipy import stats
from scipy.spatial.distance import pdist
import warnings
import pandas as pd
import sys

template_dir = "./FON_template"

class LatinHypercube:
    def __init__(self, n_samples, n_dimensions):
        """
        Initialize Latin Hypercube Sampling.
        
        Parameters:
        -----------
        n_samples : int
            Number of samples to generate
        n_dimensions : int
            Number of dimensions (variables)
        """
        self.n_samples = n_samples
        self.n_dimensions = n_dimensions
        
    def _generate_basic_lhs(self):
        """Generate basic Latin Hypercube samples without optimization."""
        # Generate the intervals
        cut = np.linspace(0, 1, self.n_samples + 1)
        
        # Create the sampling points
        samples = np.zeros((self.n_samples, self.n_dimensions))
        
        # Fill points for each dimension
        for i in range(self.n_dimensions):
            # Generate random positions within each interval
            samples[:, i] = stats.uniform(cut[:-1], np.diff(cut)).rvs(self.n_samples)
            
            # Randomly shuffle the samples
            np.random.shuffle(samples[:, i])
            
        return samples
    
    def _compute_correlation(self, samples):
        """Compute the correlation matrix of the samples."""
        return np.corrcoef(samples.T)
    
    def _objective_correlation(self, samples):
        """Calculate the sum of squares of correlation matrix elements."""
        corr = self._compute_correlation(samples)
        # Sum of squares of off-diagonal elements
        return np.sum(corr * corr) - self.n_dimensions
    
    def _maximize_min_distance(self, samples, iterations=100):
        """Maximize the minimum distance between points."""
        best_samples = samples.copy()
        best_min_dist = np.min(pdist(samples))
        
        for _ in range(iterations):
            # Randomly select dimension
            dim = np.random.randint(0, self.n_dimensions)
            
            # Randomly select two points
            i, j = np.random.choice(self.n_samples, 2, replace=False)
            
            # Swap the coordinates
            samples[i, dim], samples[j, dim] = samples[j, dim], samples[i, dim]
            
            min_dist = np.min(pdist(samples))
            if min_dist > best_min_dist:
                best_min_dist = min_dist
                best_samples = samples.copy()
            else:
                # Revert the swap
                samples[i, dim], samples[j, dim] = samples[j, dim], samples[i, dim]
                
        return best_samples
    
    def sample(self, method='basic', iterations=100, bounds=None):
        """
        Generate Latin Hypercube samples.
        
        Parameters:
        -----------
        method : str
            Sampling method: 'basic', 'corr' (correlation reduction), 
            or 'maximin' (maximize minimum distance)
        iterations : int
            Number of iterations for optimization
        bounds : array-like, shape (n_dimensions, 2)
            Lower and upper bounds for each dimension
            
        Returns:
        --------
        samples : ndarray
            The generated Latin Hypercube samples
        """
        if method not in ['basic', 'corr', 'maximin']:
            raise ValueError("Method must be 'basic', 'corr', or 'maximin'")
            
        # Generate basic LHS
        samples = self._generate_basic_lhs()
        
        # Apply optimization if requested
        if method == 'corr':
            best_samples = samples.copy()
            best_corr = self._objective_correlation(samples)
            
            for _ in range(iterations):
                # Randomly select dimension
                dim = np.random.randint(0, self.n_dimensions)
                
                # Randomly select two points
                i, j = np.random.choice(self.n_samples, 2, replace=False)
                
                # Swap the coordinates
                samples[i, dim], samples[j, dim] = samples[j, dim], samples[i, dim]
                
                corr = self._objective_correlation(samples)
                if corr < best_corr:
                    best_corr = corr
                    best_samples = samples.copy()
                else:
                    # Revert the swap
                    samples[i, dim], samples[j, dim] = samples[j, dim], samples[i, dim]
                    
            samples = best_samples
            
        elif method == 'maximin':
            samples = self._maximize_min_distance(samples, iterations)
        
        # Scale to bounds if provided
        if bounds is not None:
            bounds = np.asarray(bounds)
            if bounds.shape != (self.n_dimensions, 2):
                raise ValueError("Bounds must have shape (n_dimensions, 2)")
            samples = bounds[:, 0] + samples * (bounds[:, 1] - bounds[:, 0])
            
        return samples

def generate_lhstrainingset(seed=0, n_samples=20, n_dimensions=2):
    np.random.seed(seed)
    bounds = np.array([[-4,4] for i in range(n_dimensions)])
    dv_names = [f"x{i+1}" for i in range(n_dimensions)]

    lhs = LatinHypercube(n_samples, n_dimensions)
    samples_maximin = lhs.sample(method='maximin', bounds=bounds, iterations=1000)
    
    member_names = [f"gen=0_training={i}" for i in range(n_samples)]
    samples_df = pd.DataFrame(samples_maximin, index=member_names, columns=dv_names)
    samples_df.rename_axis("real_name", inplace=True)
    samples_df.to_csv(os.path.join(template_dir, "gp.lhs.dv_pop.csv"), index=True)

def generate_lhsstarter(seed=0, n_samples=20, n_dimensions=50):
    np.random.seed(seed+100)
    bounds = np.array([[-4,4] for i in range(n_dimensions)])
    dv_names = [f"x{i+1}" for i in range(n_dimensions)]

    lhs = LatinHypercube(n_samples, n_dimensions)
    samples_maximin = lhs.sample(method='maximin', bounds=bounds, iterations=1000)
    
    member_names = [f"gen=0_member={i}" for i in range(n_samples)]
    samples_df = pd.DataFrame(samples_maximin, index=member_names, columns=dv_names)
    samples_df.rename_axis("real_name", inplace=True)
    samples_df.to_csv(os.path.join(template_dir, "starter.dv_pop.csv"), index=True)

if __name__ == "__main__":

    seed = int(sys.argv[1])

    generate_lhstrainingset(seed)
    generate_lhsstarter(seed)
    