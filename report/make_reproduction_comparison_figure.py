import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image

def main():
    # Load images:
    # 1. Nature Ext Data Fig 1C crop
    fig1c = mpimg.imread("paper_assets/panelC_highres_crop.jpg")
    # 2. Galaxy Tutorial PI width
    gal_tut = mpimg.imread("local_reference/galaxy_smoke_pi_width.png")
    # 3. Galaxy Full Dataset 15 PI width
    gal_full = mpimg.imread("local_reference/full_dataset15_galaxy_pi_width.png")
    # 4. Galaxy Full Dataset 15 Pred vs Actual
    gal_pred = mpimg.imread("local_reference/full_dataset15_galaxy_prediction_vs_actual.png")

    fig = plt.figure(figsize=(16, 12))
    
    # Grid: top row Nature published vs Galaxy reproduction (tutorial input)
    ax1 = plt.subplot(2, 2, 1)
    ax1.imshow(fig1c)
    ax1.set_title("A. Published Extended Data Figure 1C (Nature 2026)", fontsize=13, fontweight="bold")
    ax1.axis("off")

    ax2 = plt.subplot(2, 2, 2)
    ax2.imshow(gal_tut)
    ax2.set_title("B. Galaxy UDT Reproduction (Tutorial Dataset 15 Input)", fontsize=13, fontweight="bold")
    ax2.axis("off")

    ax3 = plt.subplot(2, 2, 3)
    ax3.imshow(gal_full)
    ax3.set_title("C. Galaxy UDT Reproduction (Full 30k-Gene Dataset 15 Input)", fontsize=13, fontweight="bold")
    ax3.axis("off")

    ax4 = plt.subplot(2, 2, 4)
    ax4.imshow(gal_pred)
    ax4.set_title("D. Galaxy UDT Actual vs Predicted Acta2 (Full Dataset 15)", fontsize=13, fontweight="bold")
    ax4.axis("off")

    plt.tight_layout()
    plt.savefig("report/reproduction_comparison_panel.png", dpi=300, bbox_inches="tight")
    print("Saved report/reproduction_comparison_panel.png")

if __name__ == "__main__":
    main()
