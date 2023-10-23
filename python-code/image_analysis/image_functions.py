#!/usr/bin/env python
"""
This module stores methods required for analysis of the MP Gold CyCIF data
Authors: Nhan Huynh & Veronika Pister
"""

import math
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import cv2
import tifffile as tiff
from skimage import io, filters, exposure, morphology, measure, transform

from scipy.stats import chi2_contingency, f_oneway, norm
from scipy.stats import kstest
from scipy.stats import ttest_ind


plt.rcParams['figure.figsize'] = (10,10)


def load_image(path) -> np.ndarray:
    """
    Load image into np array variable

    Params:
        path (str) : path to the image
    Returns:
        np.ndarry
    """
    image = io.imread(path)
    blurred_image = filters.gaussian(image, sigma=1)
    return blurred_image


def stretch_image(img) -> np.ndarray:
    """
    Rescale image intensities to be visible when plotted

    Params:
        img (ndarray) : img to rescale
    Returns:
        np.ndarry
    """
    p1, p2 = np.percentile(img, (5, 95))
    img_stretched = exposure.rescale_intensity(img, in_range=(p1, p2))
    return img_stretched


def save_image(image, name):
    """
    Save np.ndarray as a .tif

    Params:
        image (ndarray) : img to save
        name (str) : name of output file
    """
    io.imsave(name, image.astype('uint8') * 255)
    
    
def display_image(image):
    """
    Plots image

    Params:
        image (ndarray) : img to plot
        
    """
    plt.imshow(image, cmap='gray')
    plt.axis('off')
    plt.show()
    
    
def load_dataset(dataset):
    """
    Loads in the MAP2, VSNL1, Taurine, and DAPI plots from each dataset

    Params:
        path (str) : path to the directory storing the dataset
    Returns:
        4 np.ndarrays
    """
    map2 = load_image(f'{dataset}/MAP2.tif')
    taurine = load_image(f'{dataset}/Taurine.tif')
    dapi = load_image(f'{dataset}/DAPI.tif')
    vsnl1 = load_image(f'{dataset}/VSNL1.tif')
    
    return map2, taurine, dapi, vsnl1


def plot_hist(img, xlim=(0,0.3)):
    """
    Plots histogram of intensity values in an img

    Params:
        img (ndarray) : img to plot
        xlim (float,float) : x lim bounds for histogram
        
    """
    with plt.rc_context({'figure.figsize': (6,6)}):
        plt.hist(img.flatten(), bins=500)
        plt.xlim(xlim)
        #plt.title("VSNL1")
        plt.xlabel("Value")
        plt.ylabel("Frequency")

        plt.grid(True)
        plt.show()
        
        
def get_mask(img, cutoff, dataset,
             structuring_element_size=2, 
             min_hole_size=20000, min_object_size=50000, 
             title="") -> np.ndarray :
    """
    Returns binary mask of positive regions based on the manual cutoff
    
    Params:
        img (np.ndarray) : the image of the marker to mask
        cutoff (float): the intensity value which qualifies a "positive" pixel
        dataset (str): title of the dataset analyzed
        structuring_element_size (int) : size of the disk to be used during closing
        min_hole_size (int) : maximum pixel size of holes in mask to be fill  
        min_object_size (int) : minimal pixel area required to keep an object
        
    Returns:
        Binary np.ndarray (mask)
    """
    structuring_element = morphology.disk(structuring_element_size)

    img_closed = morphology.binary_closing(img > cutoff, structuring_element)
    img_hole = morphology.remove_small_holes(img_closed, area_threshold=min_hole_size)
    img_obj = morphology.remove_small_objects(img_hole, min_size=min_object_size)
    
    save_image(img_obj,f'./mb{dataset}/{title}_obj.tif')
    save_image(img_hole,f'./mb{dataset}/{title}_hole.tif')
    save_image(img_closed,f'./mb{dataset}/{title}_closed.tif')
    
    return img_obj


def get_tissue_mask(img, cutoff, dataset,
                    structuring_element_size=2, 
                    min_hole_size=20000, min_object_size=50000, 
                    title=""):
    """
    Returns binary mask of actual tissue
    
    Params:
        img (np.ndarray) : the dapi image to mask
        cutoff (float): the intensity value which qualifies a "positive" pixel
        dataset (str): title of the dataset analyzed
        structuring_element_size (int) : size of the disk to be used during closing
        min_hole_size (int) : maximum pixel size of holes in mask to be fill  
        min_object_size (int) : minimal pixel area required to keep an object
        
    Returns:
        Binary np.ndarray (mask)
    """
    structuring_element = morphology.disk(structuring_element_size)

    img_closed = morphology.binary_closing(img < cutoff, structuring_element)
    img_hole = morphology.remove_small_holes(img_closed, area_threshold=min_hole_size)
    img_obj = morphology.remove_small_objects(img_hole, min_size=min_object_size)
    
    save_image(img_obj,f'./mb{dataset}/{title}_obj.tif')
    save_image(img_hole,f'./mb{dataset}/{title}_hole.tif')
    save_image(img_closed,f'./mb{dataset}/{title}_closed.tif')
    
    return np.logical_not(img_obj)


def get_taurine_by_region(tissue_mask, vsnl1_mask, map2_mask, taurine):
    """
    Returns df of location and taurine intensity for each pixel in each region.
    Takes in masks and separate image by which markers are positive in each
    region. Then quantifies the taurine dsitribution of each tissue type.
    
    Params:
        tissue_mask (np.ndarray) : binary mask array for dapi
        vsnl1_mask (np.ndarray) : binary mask array for vsnl1
        map2_mask (np.ndarray) : binary mask array for map2
        taurine (np.ndarray) : real taurine image
    """
    
    vsnl1_mask = np.logical_and(vsnl1_mask, tissue_mask)
    map2_mask = np.logical_and(map2_mask, tissue_mask)
    
    # MAKE CALLING MASKS
    map2_vsnl_obj = np.logical_and(vsnl1_mask, map2_mask)
    map2_NOTvsnl_obj = np.logical_and(np.logical_not(vsnl1_mask),map2_mask)
    vsnl_NOTmap2_obj = np.logical_and(np.logical_not(map2_mask),vsnl1_mask)
    neither_obj = np.logical_and(np.logical_not(map2_mask),np.logical_not(vsnl1_mask))
    neither_obj = np.logical_and(neither_obj, tissue_mask)
    
    # Display
    with plt.rc_context({'figure.figsize': (5,5)}):
        print("BOTH")
        display_image(map2_vsnl_obj)
        print("MAP2")
        display_image(map2_NOTvsnl_obj)
        print("VSNL1")
        display_image(vsnl_NOTmap2_obj)
        print("NEITHER")
        display_image(neither_obj)
    
    # Defining Taurine in each region
    vsnl_df = pd.DataFrame(np.log(taurine[np.logical_and(vsnl_NOTmap2_obj,taurine!=0)]),columns=["x"])
    vsnl_df["Hue"] = ["VSNL+ MAP2-" for i in range(len(vsnl_df))]

    map2_df = pd.DataFrame(np.log(taurine[np.logical_and(map2_NOTvsnl_obj,taurine!=0)]),columns=["x"])
    map2_df["Hue"] = ["MAP2+ VSNL-" for i in range(len(map2_df))]

    both_df = pd.DataFrame(np.log(taurine[np.logical_and(map2_vsnl_obj,taurine!=0)]),columns=["x"])
    both_df["Hue"] = ["MAP2+ VSNL+" for i in range(len(both_df))]

    neither_df = pd.DataFrame(np.log(taurine[np.logical_and(neither_obj,taurine!=0)]),columns=["x"])
    neither_df["Hue"] = ["MAP2- VSNL-" for i in range(len(neither_df))]

    return both_df, vsnl_df, map2_df, neither_df


def get_taurine_dist_plot(dfs, dataset):
    
    """
    Fast version of get_facet_plot (see below)
    """
    
    dfs = [df.dropna() for df in dfs]
    both_df, _, map2_df, neither_df = dfs
    
    with plt.rc_context({'figure.figsize': (14,10)}):

        fig, axs = plt.subplots(nrows=3, ncols=1,sharex=True, sharey=True)
        axs[0].hist(both_df.x, bins=500, alpha=0.7, label= "Both", color='#4285F4')
        axs[1].hist(map2_df.x, bins=500, alpha=0.7, label= "Map2+", color='#FBBC05')
        axs[2].hist(neither_df.x, bins=500, alpha=0.7, label= "Neither", color='#EA4335')
        for ax in axs:
            ax.set_xlim(-7,0)
            ax.legend()
            #ax.set_ylabel("Frequency")
        axs[-1].set_xlabel("Intensity")

        fig.text(0.04, 0.5, 'Freq.', va='center', rotation='vertical')


        #plt.legend()
        plt.suptitle("taurine distributions")

        #plt.grid(True)
        plt.savefig(f"mb{dataset}/taurine_hist.png")
        plt.show()
        
        
def get_taurine_bar_plot(dfs, dataset):
    """
    Plot the mean taurine intensity in each region
    
    Params:
        dfs: touple of four ndarrays with taurine intensities of each region
        dataset (str) : name of the dataset to be analyzed
    """
    
    dfs = [df.dropna() for df in dfs]
    both_df, _, map2_df, neither_df = dfs
    plt.clf()
    
    mean_dict = { 
    "Both (+)": np.exp(both_df.x.mean()),
    "Neither (-)": np.exp(neither_df.x.mean()),
    "MAP2+": np.exp(map2_df.x.mean())
            }

    means_df = pd.DataFrame.from_dict(mean_dict, orient="index",columns=['Dist. Mean'])
    means_df['index'] = means_df.index

    with plt.rc_context({'figure.figsize': (5,5)}):
        sns.barplot(means_df.sort_values("Dist. Mean", ascending=False), x="index", y="Dist. Mean",
                palette=['#4285F4','#FBBC05', '#EA4335'], alpha=0.8
               )
        plt.title("Mean Taurine Intensity")
        plt.savefig(f"mb{dataset}/taurine_means_bar.png")
        plt.show()


def KS_tests(dfs):
    """
    Calculate the difference between taurine distributions of different regions
    
    Params:
        dfs: touple of four ndarrays with taurine intensities of each region
    """
        
    print("KS Tests!")
    print("--------\n")
    
    dfs = [df.dropna() for df in dfs]
    both_df, _, map2_df, neither_df = dfs

    print("VSNL+/MAP+ vs MAP+/VSNL-")
    res= kstest(both_df.x,map2_df.x)
    print(f"Statistic: {round(res.statistic,2)}")
    print(f"Pval: {res.pvalue}")
    print()

    print("MAP+/VSNL- vs VSNL-/MAP-")
    res = kstest(map2_df.x,neither_df.x)
    print(f"Statistic: {round(res.statistic,2)}")
    print(f"Pval: {res.pvalue}")
    print()

    print("VSNL+/MAP+ vs VSNL-/MAP-")
    res = kstest(both_df.x,neither_df.x)
    print(f"Statistic: {round(res.statistic,2)}")
    print(f"Pval: {res.pvalue}")
    print()
    
    
def T_tests(dfs):
    
    """
    Calculate the difference in means of taurine distributions across regions
    
    Params:
        dfs: touple of four ndarrays with taurine intensities of each region
    """
    
    print("T Tests!")
    print("--------\n")
    
    dfs = [df.dropna() for df in dfs]
    both_df, _, map2_df, neither_df = dfs

    print("VSNL+/MAP+ vs MAP+/VSNL-")
    res= ttest_ind(np.exp(both_df.x),np.exp(map2_df.x))
    print(f"Statistic: {round(res.statistic,2)}")
    print(f"Pval: {res.pvalue}")
    print()

    print("MAP+/VSNL- vs VSNL-/MAP-")
    res = ttest_ind(np.exp(map2_df.x),np.exp(neither_df.x))
    print(f"Statistic: {round(res.statistic,2)}")
    print(f"Pval: {res.pvalue}")
    print()

    print("VSNL+/MAP+ vs VSNL-/MAP-")
    res = ttest_ind(np.exp(both_df.x),np.exp(neither_df.x))
    print(f"Statistic: {round(res.statistic,2)}")
    print(f"Pval: {res.pvalue}")
    print()
    
    
def get_facet_plot(both_df, map2_df, neither_df, title):
    """
    Plot the taurine distribution from each region using the same intensities
    from the x axis.
    
    Params:
        both_df (pd.DataFrame) : dataframe with the taurine intensities from 
        regions which are positive for both MAP2 and VSNL1
        map2_df (pd.DataFrame) : dataframe with the taurine intensities from 
        regions which are positive for MAP2 but not VSNL1
        neither_df (pd.DataFrame) : dataframe with the taurine intensities from
        tissue regions which are positive for neither MAP2 nor VSNL1
        title (str) : title of the distribution plot
    """
    
    concat_df = both_df.append(map2_df).append(neither_df)
    #sns.set(font_scale=2)
    
    g = sns.FacetGrid(concat_df, #the dataframe to pull from
                  row="Hue", #define the column for each subplot row to be differentiated by
                  hue="Hue", #define the column for each subplot color to be differentiated by
                  aspect=8, #aspect * height = width
                  height=3, #height of each subplot
                  palette=['#4285F4','#EA4335','#FBBC05'] #google colors
                 )
    
    g.map(sns.kdeplot, "x", fill=True, alpha=1, lw=1.5, bw_method=0.2)
    g.map(sns.kdeplot, "x", lw=4, bw_method=0.2)
    
    g.fig.subplots_adjust(top=0.875, bottom=0.1) 
    g.fig.suptitle(title + " Taurine Distributions by Region", fontsize=30)
    g.axes[2,0].set_xlabel('pixel intensity (log scaled)', fontsize = 20)
    
    g.axes[2,0].set_ylabel('', fontsize = 25)
    g.axes[0,0].set_ylabel('', fontsize = 25)
    g.axes[1,0].set_ylabel('Density', fontsize = 20)

    g.axes[0,0].set_title('MAP2+ VSNL1+', fontsize = 22.5, y=0.9)
    g.axes[1,0].set_title('MAP2+ VSNL1-', fontsize = 22.5, y=0.9)
    g.axes[2,0].set_title('MAP2- VSNL1-', fontsize = 22.5, y=0.9)
    
    g.tight_layout()
    g.savefig(title.replace("-","_") + "_facet.svg")
    
    return g