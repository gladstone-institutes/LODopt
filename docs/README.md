# LODopt Website

This directory contains the GitHub Pages website for the LODopt R package.

## Files

- `index.html` - Main landing page with overview, installation, and quick start guide
- `intro.html` - Tutorial page (rendered from vignettes/intro.Rmd)
- `documentation.html` - Comprehensive documentation with function details and examples
- `README.md` - This file

## GitHub Pages Setup

To enable GitHub Pages for this repository:

1. Go to your repository settings on GitHub
2. Scroll down to the "Pages" section
3. Under "Source", select "Deploy from a branch"
4. Choose the branch (usually `main` or `master`)
5. Set the folder to `/docs`
6. Click "Save"

The website will be available at `https://[username].github.io/LODopt/`

## Website Structure

The website follows a similar structure to other Gladstone Institutes packages like clustOpt:

- **Home page** (`index.html`): Overview, installation instructions, and quick start
- **Tutorial** (`intro.html`): Step-by-step guide with examples
- **Documentation** (`documentation.html`): Comprehensive function documentation and troubleshooting

## Styling

The website uses a clean, modern design with:
- Responsive layout that works on desktop and mobile
- Consistent color scheme matching Gladstone Institutes branding
- Clear navigation between pages
- Syntax-highlighted code blocks
- Professional typography

## Maintenance

To update the website:
1. Edit the HTML files in this directory
2. Commit and push changes to GitHub
3. GitHub Pages will automatically rebuild the site

For the tutorial page (`intro.html`), update the source R Markdown file in `vignettes/intro.Rmd` and re-render it.
