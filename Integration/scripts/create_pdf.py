from PIL import Image
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
import math
from pypdf import PdfWriter
import argparse
import glob
from PIL import Image, ImageDraw, ImageFont
import matplotlib.pyplot as plt
import os
import cv2
import PIL
import cairosvg
import tempfile
from xml.etree import ElementTree as ET
import math

from reportlab.pdfgen import canvas
from reportlab.lib.units import mm
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF
from reportlab.lib.pagesizes import A4
import cairosvg

def convert_svg_to_pdf(svg_path, pdf_path):
    """
    Converts an SVG file to a PDF using CairoSVG.
    """
    cairosvg.svg2pdf(url=svg_path, write_to=pdf_path)
    print(f"Converted {svg_path} to {pdf_path}")
    
    
PAPERLESS_OCR_MAX_IMAGE_PIXELS=40000000000
Image.MAX_IMAGE_PIXELS = None

def textsize(text, font):
    im = Image.new(mode="P", size=(0, 0))
    draw = ImageDraw.Draw(im)
    _, _, width, height = draw.textbbox((0, 0), text=text, font=font)
    return width, height



def create_pdf_from_single_svg(svg_path, output_pdf_path, page_size=A4, margin=50):
    """
    Renders a single SVG onto the center of a PDF page, scaled to fit while preserving aspect ratio.
    """
    page_width, page_height = page_size
    c = canvas.Canvas(output_pdf_path, pagesize=page_size)

    drawing = svg2rlg(svg_path)
    if not drawing or getattr(drawing, "width", 0) == 0 or getattr(drawing, "height", 0) == 0:
        raise ValueError(f"{svg_path} has invalid or missing dimensions.")

    # Available space for the drawing
    available_width = page_width - 2 * margin
    available_height = page_height - 2 * margin

    # Compute scaling factor (fit within available space)
    scale_x = available_width / drawing.width
    scale_y = available_height / drawing.height
    scale = min(scale_x, scale_y)

    # Compute final size
    render_width = drawing.width * scale
    render_height = drawing.height * scale

    # Centered position
    x = (page_width - render_width) / 2
    y = (page_height - render_height) / 2

    # Draw
    c.saveState()
    c.translate(x, y)
    c.scale(scale, scale)
    renderPDF.draw(drawing, c, 0, 0)
    c.restoreState()

    c.save()
    print(f"Saved single-SVG PDF to {output_pdf_path}")
    

def create_pdf_with_vector_svg_grid_and_captions(svg_paths, captions, header_text, output_pdf_path,
                                                 n_cols=1, page_size=A4, margin=20, spacing=20,
                                                 caption_height=30, added=False, fallback_size=(1500, 1500)):
    """
    Create a PDF with vector-based SVGs arranged in a grid with captions and a header.
    Scales each SVG uniformly to fit a defined grid cell to avoid overlap.
    """
    if len(svg_paths) != len(captions):
        raise ValueError("Number of captions must match number of SVGs.")

    page_width, page_height = page_size
    c = canvas.Canvas(output_pdf_path, pagesize=page_size)

    # Draw header
    c.setFont("Helvetica-Bold", 16)
    c.drawCentredString(page_width / 2, page_height - margin, header_text)

    start = 0
    if added:
        svg_subset = svg_paths[:2]
        start = 2
    else:
        svg_subset = []
    svg_subset.extend(svg_paths[start:])

    n_images = len(svg_subset)
    n_rows = math.ceil(n_images / n_cols)

    grid_width = page_width - 2 * margin
    grid_height = page_height - 2 * margin - 40

    cell_width = (grid_width - (n_cols - 1) * spacing) / n_cols
    cell_height = (grid_height - (n_rows - 1) * spacing) / n_rows
    image_height_max = cell_height - caption_height

    start_x = margin
    start_x = 0
    start_y = page_height - margin - 40
    
    print("svg_subset")
    print(svg_subset)
    print(start_y)
    print(start_x)

    for idx, (svg_path, caption) in enumerate(zip(svg_subset, captions)):
        row = idx // n_cols
        col = idx % n_cols

        x = start_x + col * (cell_width + spacing)
        y = start_y - row * (cell_height + spacing) - image_height_max - caption_height

        drawing = svg2rlg(svg_path)

        # Handle missing size info
        if getattr(drawing, 'width', None) is None or getattr(drawing, 'height', None) is None:
            print(f"Warning: {svg_path} missing size info. Forcing fallback size.")
            drawing.width, drawing.height = fallback_size

        #drawing.width, drawing.height = fallback_size
        
        drawing_width = drawing.width
        drawing_height = drawing.height
        
        print("SIZES")
        print(drawing_width)
        print(drawing_height)

        # Uniform scaling
        print("SCALING")
        scale_x = cell_width / drawing_width
        scale_y = image_height_max / drawing_height
        scale = min(scale_x, scale_y) * 0.3
        print(scale_x)
        print(scale_y)
        print(scale)

        scaled_width = drawing_width * scale
        scaled_height = drawing_height * scale
        
        print("cell_width")
        print(cell_width)
        print(x)

        draw_x = x + (cell_width - scaled_width) / 2
        draw_y = y + caption_height + (image_height_max - scaled_height) / 2

        c.saveState()
        c.translate(draw_x, draw_y)
        c.scale(scale, scale)
        try:
            renderPDF.draw(drawing, c, 0, 0)
        except Exception as e:
            print(f"Failed to render {svg_path}: {e}")
        c.restoreState()

        # Draw caption
        c.setFont("Helvetica", 10)
        c.drawCentredString(x + cell_width / 2, y, caption)

    c.save()
    print(f"PDF with vector SVGs saved to {output_pdf_path}")

    
def create_image_grid_with_captions(png_paths, captions, header_text, output_png_path,
                                    n_cols=3, image_size=(800, 1000), margin=40,
                                    spacing=10, caption_height=20, header_height=40, added=False):
    """
    Combine images into a single PNG arranged in a grid with captions and a header.

    Args:
        png_paths (list of str): List of PNG image file paths.
        captions (list of str): Captions for each image.
        header_text (str): Header text for the image.
        output_png_path (str): Output PNG file path.
        n_cols (int): Number of columns in the grid.
        image_size (tuple): Size of the final image (width, height).
        margin (int): Margin around content.
        spacing (int): Spacing between images.
        caption_height (int): Height of the caption area below each image.
        header_height (int): Height reserved for the header text.
        added (bool): Whether to treat first two elements differently.
    """
    
    header_height = 80  # Space for top header
    image_title_height = 40  # Space above each image
    padding = 20  # Between elements
    
    
    if len(png_paths) != len(captions):
        raise ValueError("Number of captions must match number of images.")

    total_width, total_height = image_size

    start = 0
    if added:
        images = [Image.open(p).convert("RGB") for p in png_paths[:2]]
        start = 2
    else:
        images = []

    images.extend([Image.open(p).convert("RGB") for p in png_paths[start:]])
    n_images = len(images)
    n_rows = math.ceil(n_images / n_cols)

    grid_width = total_width - 2 * margin
    grid_height = total_height - 2 * margin - header_height

    cell_width = (grid_width - (n_cols - 1) * spacing) / n_cols
    cell_height = (grid_height - (n_rows - 1) * spacing) / n_rows

    image_max_height = cell_height - caption_height

    # Create the final image
    result = Image.new("RGB", (total_width, total_height), "white")
    draw = ImageDraw.Draw(result)

    try:
        font = ImageFont.truetype("arial.ttf", size=32)  # image title
        header_font = ImageFont.truetype("arial.ttf", size=56)  # page header
    except:
        print("⚠️ Could not load Arial font. Using default font.")
        font = ImageFont.load_default(size=32)
        header_font = ImageFont.load_default(size=56)

    # Draw header
    header_y = margin
    header_x = total_width // 2
    draw.text((header_x, header_y), header_text, font=header_font, anchor="mm")

    start_x = margin
    start_y = margin + header_height

    for idx, (img, caption) in enumerate(zip(images, captions)):
        row = idx // n_cols
        col = idx % n_cols

        x = start_x + col * (cell_width + spacing)
        y = start_y + row * (cell_height + spacing)

        scale = min(cell_width / img.width, image_max_height / img.height)
        new_width = max(1, int(img.width * scale))
        new_height = max(1, int(img.height * scale))

        resized_img = img.resize((new_width, new_height), Image.LANCZOS)

        # Position image centered in its cell
        img_x = int(x + (cell_width - new_width) / 2)
        img_y = int(y)

        result.paste(resized_img, (img_x, img_y))

        # Draw caption
        caption_x = x + cell_width / 2
        caption_y = y + new_height + 2
        draw.text((caption_x, caption_y), caption, font=font, fill="black", anchor="mm")

    result.save(output_png_path)
    print(f"Grid image with captions saved to {output_png_path}")
    
    
def create_pdf_with_image_grid_and_captions(png_paths, captions, header_text, output_pdf_path, n_cols=3, page_size=A4, margin=40, spacing=10, caption_height=20, added = False):
    """
    Create a PDF with images arranged in a grid, a header, and captions under each image.

    Args:
        png_paths (list of str): List of PNG image file paths.
        captions (list of str): Captions for each image.
        header_text (str): Header text for the page.
        output_pdf_path (str): Output PDF file path.
        n_cols (int): Number of columns in the grid.
        page_size (tuple): PDF page size (e.g., A4).
        margin (int): Margin around content.
        spacing (int): Spacing between images.
        caption_height (int): Vertical space reserved for each caption.
    """
    if len(png_paths) != len(captions):
        raise ValueError("Number of captions must match number of images.")

    page_width, page_height = page_size
    c = canvas.Canvas(output_pdf_path, pagesize=page_size)

    # Header
    c.setFont("Helvetica-Bold", 16)
    c.drawCentredString(page_width / 2, page_height - margin, header_text)

    # Load and prepare images
    start = 0
    if added == True: 
        images = [png_paths[0], png_paths[1]]
        start = 2
    else: images = []
    images.extend([Image.open(p).convert("RGB") for p in png_paths[start:]])
    
    print(images)
    
    n_images = len(images)
    n_rows = math.ceil(n_images / n_cols)

    # Grid dimensions
    grid_width = page_width - 2 * margin
    grid_height = page_height - 2 * margin - 40  # Leave room for header

    cell_width = (grid_width - (n_cols - 1) * spacing) / n_cols
    cell_height = (grid_height - (n_rows - 1) * spacing) / n_rows

    image_height_max = cell_height - caption_height

    start_x = margin
    start_y = page_height - margin - 40  # below header

    for idx, (img, caption) in enumerate(zip(images, captions)):
        row = idx // n_cols
        col = idx % n_cols

        x = start_x + col * (cell_width + spacing)
        y = start_y - row * (cell_height + spacing) - image_height_max - caption_height

        # Resize image to fit inside (cell_width, image_height_max)
        print(img)
        #img.thumbnail((cell_width, image_height_max), Image.LANCZOS)
        
        # Compute scale factor (preserve aspect ratio)
        scale = min(cell_width / img.width, image_height_max / img.height)
        new_width = max(1, int(img.width * scale))
        new_height = max(1, int(img.height * scale))

        # Resize using a copy to avoid mutating original
        resized_img = img.resize((new_width, new_height), resample=Image.LANCZOS)

        # Ensure the image is RGB (some PNGs might be RGBA or P)
        if resized_img.mode not in ("RGB", "L"):
            resized_img = resized_img.convert("RGB")

        # Center image within its cell
        img_x = x + (cell_width - new_width) / 2
        img_y = y + caption_height

        # Draw resized image
        c.drawImage(ImageReader(resized_img), img_x, img_y, width=new_width, height=new_height)
        # Draw caption centered under image
        c.setFont("Helvetica", 10)
        c.drawCentredString(x + cell_width / 2, y, caption)

    c.save()
    print(f"PDF with captions saved to {output_pdf_path}")
    
    
    
def create_pdf_with_6_images(image_paths, image_titles, page_header, output_pdf_path):
    """
    Combines 6 PNG images into a single-page PDF with a large page header and image-specific headers.

    Layout: 2 images per row, 3 rows (no centering).

    Parameters:
    - image_paths: List of 6 PNG file paths.
    - image_titles: List of 6 strings, one title per image.
    - page_header: String for the page-level header.
    - output_pdf_path: Output path for the PDF file.
    """
    assert len(image_paths) == 6, "Exactly 6 image paths are required."
    assert len(image_titles) == 6, "Exactly 6 image titles are required."

    # Load and resize images to a target width (keeps aspect ratio)
    target_width = 2000  # Adjust if needed
    images = [Image.open(p).convert("RGB") for p in image_paths]
    resized_images = [
        img.resize((target_width, int(img.height * target_width / img.width)))
        for img in images
    ]
    max_height = max(img.height for img in resized_images)

    # Fonts (try Arial, then DejaVu, then default)
    def _load_font(pref, size):
        try:
            return ImageFont.truetype(pref, size=size)
        except Exception:
            return None

    font = (
        _load_font("arial.ttf", 32)
        or _load_font("DejaVuSans.ttf", 32)
        or ImageFont.load_default()
    )
    header_font = (
        _load_font("arial.ttf", 56)
        or _load_font("DejaVuSans.ttf", 56)
        or ImageFont.load_default()
    )

    # Layout
    padding = 40
    image_title_height = 80
    header_height = 120
    num_cols = 2

    total_width = num_cols * target_width + (num_cols + 1) * padding
    # header area + top padding + 3 rows * (title + image + padding) + bottom padding
    total_height = header_height + padding + 3 * (image_title_height + max_height + padding) + padding

    # Create canvas
    canvas = Image.new("RGB", (total_width, total_height), "white")
    draw = ImageDraw.Draw(canvas)

    # Draw page header (centered)
    header_bbox = draw.textbbox((0, 0), page_header, font=header_font)
    header_w = header_bbox[2] - header_bbox[0]
    draw.text(((total_width - header_w) // 2, padding // 2), page_header, fill="black", font=header_font)

    # Paste images with titles (2 per row, 3 rows)
    for idx, (img, title) in enumerate(zip(resized_images, image_titles)):
        row = idx // 2            # 0, 1, 2
        col = idx % 2             # 0 or 1
        x = padding + col * (target_width + padding)
        y = header_height + padding + row * (max_height + image_title_height + padding)

        # Title (centered above image)
        title_bbox = draw.textbbox((0, 0), title, font=font)
        title_w = title_bbox[2] - title_bbox[0]
        draw.text((x + (target_width - title_w) // 2, y), title, fill="black", font=font)

        # Image below title
        canvas.paste(img, (x, y + image_title_height))

    # Save to PDF
    canvas.save(output_pdf_path, "PDF", resolution=100.0)
    
    
def create_fig(sequence, window, coverage, name_):
    

    
    x_list = sequence

    
    color = []
    x_list_ = []
    for enum, x_ in enumerate(x_list):
            
        if x_ == "C":
            color.append("blue")
        elif x_ == "G":
            color.append("yellow")
        elif x_ == "T":
            color.append("grey")
        elif x_ == "A":
            color.append("red")
        x_list_.append(x_)
                
    x_list = x_list_
    

    tiltes = ["1lbeta-treated","untreated","RNA seq"]
    
    fig, axes = plt.subplots(1, 1, sharey=False, figsize = (30, 5))
    

            
    axes.bar(range(0, len(x_list)), window, color = color)
    axes.set_xticks(range(0, len(x_list)))
    axes.set_xticklabels(x_list, fontsize = 6)


    axes.set_title(name_)           

    axes.set_ylabel('DMS signal')           
    axes.set_ylim(bottom=0, top=0.35)        
        
    ax_twin = axes.twinx()
    ax_twin.plot(coverage[enum], color='green', linewidth=5)
    ax_twin.grid(None)
    ax_twin.set_ylabel('coverage')
    ax_twin.set_ylim(0, 5)
    ax_twin.tick_params(axis = "y", labelsize = 6)
    axes.tick_params(axis = "y", labelsize = 6)
      
    labels = ["A", "C", "G", "T", "coverage"]
    colors = ["red", "blue", "yellow", "grey", "green"]
          
    handles = [plt.Rectangle((0,0),1,1, color=colors[en]) for en,label in enumerate(labels)]
    axes.legend(handles, labels)
        
        
    fig.suptitle(name_)
    
        
    #plt.savefig(save_path + str(name_) +"_" +str(strand) + ".svg")
    #plt.close()
        
    pil_img = save_plot_and_get(fig)
    return pil_img
    
def save_plot_and_get(fig):
    fig.savefig("test.jpg")
    img = cv2.imread("test.jpg")
    os.remove("test.jpg")
    return PIL.Image.fromarray(img)
    
    
def create_text_image(text, width=800, height=400, font_size=40, font_path=None):
    """
    Create a white PNG image with centered text.

    Args:
        text (str): Text to write on the image.
        width (int): Width of the image in pixels.
        height (int): Height of the image in pixels.
        font_size (int): Font size for the text.
        font_path (str): Optional path to a .ttf font file.
    """
    # Create white background
    img = Image.new('RGB', (width, height), color='white')
    draw = ImageDraw.Draw(img)

    # Load font (default or from path)
    if font_path:
        font = ImageFont.truetype(font_path, font_size)
    else:
        font = ImageFont.load_default()

    # Calculate text size and position
    text_lines = text.split('\n')
    y_offset = (height - font_size * len(text_lines)) // 2

    for i, line in enumerate(text_lines):
        text_width, text_height = textsize(line, font=font)
        x = (width - text_width) // 2
        y = y_offset + i * font_size
        draw.text((x, y), line, fill='black', font=font)

    return img
    
def merge_pdfs(pdf_paths, output_path):
    """
    Merge multiple PDF files into one.

    Args:
        pdf_paths (list of str): List of PDF file paths to merge.
        output_path (str): Output path for the merged PDF.
    """
    print("MERGING PDFS")
    print(pdf_paths)
    merger = PdfWriter()
    for pdf in pdf_paths:
        merger.append(pdf)
    merger.write(output_path)
    merger.close()
    print(f"Merged PDF saved to {output_path}")



def create_pdf_with_vector_svg_and_text(svg_paths, text_blocks, header_text, output_pdf_path,
                                        page_size=(2000,2000), margin=40, spacing=40,
                                        svg_width_ratio=0.5, row_height=150, fallback_size=(200, 200)):
    """
    Create a PDF with SVGs on the left and corresponding text on the right.

    Args:
        svg_paths (list of str): List of SVG file paths.
        text_blocks (list of str): Text blocks to go beside each SVG.
        header_text (str): Title text at the top.
        output_pdf_path (str): Output path for the PDF.
        page_size (tuple): PDF page size.
        margin (int): Page margin in points.
        spacing (int): Spacing between rows.
        svg_width_ratio (float): Fraction of row width occupied by the SVG.
        row_height (int): Height of each row (SVG + text).
        fallback_size (tuple): Fallback size for broken SVGs.
    """
    if len(svg_paths) != len(text_blocks):
        raise ValueError("Each SVG must have a corresponding text block.")
        
    print(page_size)

    page_width, page_height = page_size
    c = canvas.Canvas(output_pdf_path, pagesize=page_size)

    # Header
    c.setFont("Helvetica-Bold", 16)
    c.drawCentredString(page_width / 2, page_height - margin, header_text)

    max_rows_per_page = int((page_height - 2 * margin - 40) / (row_height + spacing))
    current_y = page_height - margin - 40

    for idx, (svg_path, text) in enumerate(zip(svg_paths, text_blocks)):
        if idx > 0 and idx % max_rows_per_page == 0:
            c.showPage()
            c.setFont("Helvetica-Bold", 16)
            c.drawCentredString(page_width / 2, page_height - margin, header_text)
            current_y = page_height - margin - 40

        drawing = svg2rlg(svg_path)

        drawing_width = getattr(drawing, 'width', fallback_size[0])
        drawing_height = getattr(drawing, 'height', fallback_size[1])

        # Layout areas
        row_top = current_y
        row_bottom = current_y - row_height
        svg_area_width = svg_width_ratio * (page_width - 2 * margin)
        text_area_x = margin + svg_area_width + spacing
        text_area_width = page_width - margin - text_area_x

        # Scale SVG to fit in left box
        scale = min(svg_area_width / drawing_width, row_height / drawing_height)
        svg_x = margin + (svg_area_width - drawing_width * scale) / 2
        svg_y = row_bottom + (row_height - drawing_height * scale) / 2

        c.saveState()
        c.translate(svg_x, svg_y)
        c.scale(scale, scale)
        renderPDF.draw(drawing, c, 0, 0)
        c.restoreState()

        # Draw corresponding text block to the right
        c.setFont("Helvetica", 5)
        text_lines = text.split('\n')
        text_y = row_top - 10


        for line in text_lines:
            print(line)
            words = line.split('**')
            is_bold = False
            cursor_x = text_area_x
            for segment in words:
                if is_bold:
                    c.setFont("Helvetica-Bold", 5)
                else:
                    c.setFont("Helvetica", 5)
                c.drawString(cursor_x, text_y, segment)
                text_width = c.stringWidth(segment, "Helvetica-Bold" if is_bold else "Helvetica", 5)
                cursor_x += text_width
                is_bold = not is_bold
            text_y -= 10

        current_y -= row_height + spacing

    c.save()
    print(f"PDF with SVGs and text saved to {output_pdf_path}")
    
    
def fix_svg_dimensions(svg_path, output_path=None, width="500", height="500"):
    tree = ET.parse(svg_path)
    root = tree.getroot()
    
    root.set("width", width)
    root.set("height", height)
    root.attrib.pop("viewBox", None)  # optional: remove viewBox if needed

    if output_path is None:
        output_path = svg_path
    tree.write(output_path)





def create_pdf_with_5_images(image_paths, image_titles, page_header, output_pdf_path):
    """
    Combines 5 PNG images into a single-page PDF with a large page header and image-specific headers.

    Layout: 2 images in first two rows, 1 centered image in third row.

    Parameters:
    - image_paths: List of 5 PNG file paths.
    - image_titles: List of 5 strings, one title per image.
    - page_header: String for the page-level header.
    - output_pdf_path: Output path for the PDF file.
    """
    assert len(image_paths) == 5, "Exactly 5 image paths are required."
    assert len(image_titles) == 5, "Exactly 5 image titles are required."

    # Load and resize images to target width
    target_width = 2000  # Adjust if needed
    images = [Image.open(p).convert("RGB") for p in image_paths]
    resized_images = [
        img.resize((target_width, int(img.height * target_width / img.width)))
        for img in images
    ]
    max_height = max(img.height for img in resized_images)

    # Fonts
    try:
        font = ImageFont.truetype("arial.ttf", size=32)  # image title
        header_font = ImageFont.truetype("arial.ttf", size=56)  # page header
    except:
        print("⚠️ Could not load Arial font. Using default font.")
        font = ImageFont.load_default(size=70)
        header_font = ImageFont.load_default(size=100)

    # Layout
    padding = 40
    image_title_height = 80
    header_height = 120
    num_cols = 2
    total_width = num_cols * target_width + (num_cols + 1) * padding
    total_height = (
        header_height +
        2 * (max_height + image_title_height + padding) +  # rows 1 & 2
        (max_height + image_title_height + padding) +      # row 3 (1 image)
        padding  # bottom margin
    )

    # Create canvas
    canvas = Image.new("RGB", (total_width, total_height), "white")
    draw = ImageDraw.Draw(canvas)

    # Draw page header
    header_bbox = draw.textbbox((0, 0), page_header, font=header_font)
    header_w = header_bbox[2] - header_bbox[0]
    draw.text(((total_width - header_w) // 2, padding // 2), page_header, fill="black", font=header_font)

    # Paste images with titles
    for idx, (img, title) in enumerate(zip(resized_images, image_titles)):
        if idx < 4:
            row = idx // 2
            col = idx % 2
            x = padding + col * (target_width + padding)
        else:
            # Last image centered
            x = (total_width - target_width) // 2
            row = 2
        y = (
            header_height + padding +
            row * (max_height + image_title_height + padding)
        )

        # Draw image title
        title_bbox = draw.textbbox((0, 0), title, font=font)
        title_w = title_bbox[2] - title_bbox[0]
        draw.text((x + (target_width - title_w) // 2, y), title, fill="black", font=font)

        # Paste image below title
        canvas.paste(img, (x, y + image_title_height))

    # Save to PDF
    canvas.save(output_pdf_path, "PDF", resolution=100.0)
    
    
if __name__ == "__main__":


            
            
    cmdline_parser = argparse.ArgumentParser('retrieve RNA seq varna structure')

    cmdline_parser.add_argument('-f', '--output_folder',
                                default="",
                                help='output folder',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-c', '--dms_analysis_file',
                                default="",
                                help='dms file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-d', '--output_pdf',
                                default="",
                                help='output pdf file',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-e', '--pdf_folder',
                                default="",
                                help='pdf output folder',
                                required = True,
                                type=str)
    cmdline_parser.add_argument('-a', '--run_type',
                                default="3utr",
                                help='run_type',
                                required = True,
                                type=str)
                                
                                
    args, unknowns = cmdline_parser.parse_known_args()

    PAPERLESS_OCR_MAX_IMAGE_PIXELS=40000000000
    dms_file_open = open(args.dms_analysis_file).readlines()
    dms_dict = {}
    
    for enum, line_ in enumerate(dms_file_open):
        if line_[0] == ">":
            dms_dict[line_[1:-1]] = [dms_file_open[enum+1], dms_file_open[enum+2], dms_file_open[enum+3], dms_file_open[enum+4]]
    
    
    folders = glob.glob(args.output_folder + "/*")
    
    for folder in folders:
    
        id_ = folder.split("/")[-1]
        if id_.split(".")[-1] == "fa": continue
        header = open(args.output_folder + id_ +f"/complete_info_{args.run_type}.txt").readlines()[0]
        
        varna_files = [args.output_folder +   id_ + "/plain_RNAfold_varna.png", args.output_folder +   id_ + "/dms_RNAfold_varna.png", args.output_folder +   id_ + "/eclip_RNAfold_varna.png" , args.output_folder +   id_ + "/oops_single_RNAfold_varna.png", args.output_folder +   id_ + "/ml_dms_RNAfold_varna.png" , args.output_folder +   id_ + "/oops_single_dms_RNAfold_varna.png"]
        
        captions_utr = ["plain", "dms", "protein", "protein+dms", "dms+ml", "fully integrated"]

        try:
            create_pdf_with_6_images(varna_files, captions_utr, header, args.pdf_folder + "/" + id_  + "_utr" + ".pdf")
        except:
            max_width = max(img.width for img in varna_files)
            resized_images = [img.resize((max_width, int(img.height * max_width / img.width))) for img in varna_files]
            create_pdf_with_6_images(resized_images, captions_utr, header, args.pdf_folder + "/" + id_  + "_utr" + ".pdf")
    
    merge_pdfs(glob.glob(args.pdf_folder + "/*.pdf"), args.output_pdf)


