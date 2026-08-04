# Content Update Guide: Mesele Tilahun Belete Portfolio

Your new website is designed to be **data-driven**. This means you don't need to touch any HTML code to add new research, publications, or software. You only need to edit the `data/data.json` file.

## 1. Locating the Data File
Inside your website folder, navigate to:
`data/data.json`

Open this file in any text editor (VS Code, Notepad++, or even GitHub's online editor).

---

## 2. Adding a New Publication
To add a new publication, find the `"publications": [` section and add a new block at the **top** of the list (so it appears first on your site).

### Example Block:
```json
{
  "title": "Your New Paper Title",
  "journal": "Nature Biotechnology",
  "year": "2026",
  "doi": "10.1038/s41587-026-xxxx-x",
  "authors": "Belete, M.T., et al.",
  "type": "First Author"
}
```

**Steps:**
1. Copy an existing publication block (everything between `{` and `}`).
2. Paste it at the top of the list.
3. Update the fields. **Note:** Ensure there is a comma `,` between blocks, but **no comma** after the last block in the list.

---

## 3. Updating Software Tools
Find the `"software": [` section. Each software item is a "card" on your Software page.

### Example Block:
```json
{
  "name": "New Tool Name",
  "description": "A short one-sentence description of what the tool does.",
  "features": ["Feature 1", "Feature 2", "Feature 3"],
  "github_url": "https://github.com/mesele1621/new-tool",
  "screenshot": "images/software/new-tool.jpg",
  "status": "Beta" 
}
```
*   **Status:** You can add `"status": "Coming Soon"` or `"status": "New"` to show a badge on the card.
*   **Features:** This is a list of tags that will appear on the card.

---

## 4. Updating Research Projects
Find the `"projects": [` section. These appear on your Projects page with progress bars.

### Example Block:
```json
{
  "title": "Novel Virus Discovery in Tropical Crops",
  "description": "Investigating emerging pathogens in sub-Saharan Africa.",
  "institution": "APQA / International Collaborators",
  "period": "2026 – Present",
  "status": "Ongoing",
  "impact": "Identification of 3 novel species"
}
```
*   **Status:** Use `"Ongoing"` or `"Completed"`. 
*   **Progress Bar:** The website automatically sets the progress bar to 75% for "Ongoing" and 100% for "Completed".

---

## 5. Updating Profile Information
At the very top of the file, you can update your general details:
*   **`roles`**: These are the words that "type" out in the hero section on the home page.
*   **`stats`**: Update your citation count, h-index, or total publications here.

---

## Important Tips for Success
1.  **JSON Syntax:** JSON is strict. Always ensure your text is inside double quotes `" "`.
2.  **Commas:** Every item in a list must have a comma after it, **except the last one**.
3.  **Images:** If you add a new software screenshot, place the image in `images/software/` and make sure the path in the JSON matches exactly.
4.  **Refresh:** After saving the file and uploading it to GitHub, refresh your browser. You may need to wait 1-2 minutes for GitHub Pages to update.

> **Pro Tip:** If the website stops loading after an update, you likely have a missing comma or a stray quote in your JSON. Use a [JSON Validator](https://jsonlint.com/) to find the error quickly.
