## ppg-import

Workflow using the R [targets](https://github.com/ropensci/targets) package to
import data from the
[World Ferns](https://www.worldplants.de/world-ferns/ferns-and-lycophytes-list)
list by Michael Hassler and convert it to
[Darwin Core](https://dwc.tdwg.org/terms/#taxon) format.

The pipeline is triggered automatically when a collaborator uploads a new
`WorldFerns_ver_*.csv` file to a shared private Google Drive folder. A Google
Apps Script watches the folder and fires a GitHub repository dispatch event;
GitHub Actions then downloads the file using a Google service account, runs
the pipeline, and uploads `wf_dwc.csv` as a workflow artifact.

---

## Automated pipeline setup

### 1. Google Drive folder

1. Go to [Google Drive](https://drive.google.com) and click
   **New → New folder**. Name it something like `World Ferns uploads`.
2. Right-click the folder and choose **Share**. Add your collaborator's
   email with **Editor** access so they can upload files to it.
3. Ask your collaborator to name their files following the pattern
   `WorldFerns_ver_YY-MM.csv` (e.g. `WorldFerns_ver_26-05.csv`) so the
   pipeline can detect them automatically. Files with other names in the
   same folder are ignored.

### 2. Create a GitHub Personal Access Token (PAT)

The Apps Script bridge needs a token to tell GitHub to start a workflow run.

1. Go to [github.com/settings/tokens](https://github.com/settings/tokens)
   (you must be logged in to GitHub as an owner or collaborator of
   `pteridogroup/ppg-import`).
2. Click **Generate new token → Generate new token (classic)**.
3. Give it a descriptive name, e.g. `ppg-import Apps Script trigger`.
4. Set the **Expiration** to something reasonable (90 days or 1 year).
   You will need to rotate it before it expires.
5. Under **Select scopes**, tick **`repo`** (the full repo scope is needed
   to fire `repository_dispatch` events on a private repo). If the repo is
   public, `public_repo` is sufficient.
6. Click **Generate token**.
7. **Copy the token immediately** — GitHub will only show it once. Paste
   it somewhere safe (e.g. a password manager) until you add it to Apps
   Script in step 3d.

### 3. Google Apps Script bridge

Google Apps Script is a free JavaScript scripting environment built into Google
Drive. You do not need to install anything — it runs entirely in the browser.

#### 3a. Find your Drive folder ID

1. Open [Google Drive](https://drive.google.com) and navigate to the folder
   where your collaborator will upload files.
2. Look at the browser URL bar. It will look like:
   ```
   https://drive.google.com/drive/folders/1A2B3C4D5E6F7G8H9I0J
   ```
3. Copy the long string of letters and numbers at the end — that is your
   **folder ID**. Keep it handy for step 3d.

#### 3b. Create an Apps Script project

1. With the Drive folder open, click **New → More → Google Apps Script**.
   (If you don't see Apps Script, click **Connect more apps**, search for
   "Apps Script", and install it. Then repeat.)
2. A new browser tab opens with the Apps Script editor.
3. At the top, click the project name ("Untitled project") and rename it
   to something like **ppg-import trigger**.

#### 3c. Paste the script code

1. In the editor you will see a file called `Code.gs` with a default empty
   function. **Select all the existing text and delete it.**
2. Paste the entire block below:

```javascript
// ── Configuration (values are stored in Script Properties, not here) ── //
const GITHUB_TOKEN    = PropertiesService
  .getScriptProperties().getProperty("GITHUB_TOKEN");
const GITHUB_REPO     = "pteridogroup/ppg-import";
const DRIVE_FOLDER_ID = PropertiesService
  .getScriptProperties().getProperty("DRIVE_FOLDER_ID");
// ───────────────────────────────────────────────────────────────────── //

/**
 * Run this function ONCE from the editor to install the Drive onChange
 * trigger. After that, onDriveChange() fires automatically whenever any
 * file is added to your Drive (the function ignores irrelevant files).
 * Re-running installDriveTrigger() is safe — it removes any existing
 * copy first to avoid duplicates.
 */
function installDriveTrigger() {
  // Remove any pre-existing triggers for onDriveChange.
  ScriptApp.getProjectTriggers().forEach(function(t) {
    if (t.getHandlerFunction() === "onDriveChange") {
      ScriptApp.deleteTrigger(t);
    }
  });
  ScriptApp.newTrigger("onDriveChange")
    .forDrive()
    .onChange()
    .create();
  Logger.log("Drive onChange trigger installed.");
}

/**
 * Called automatically by the Drive onChange trigger whenever a file
 * is created, modified, or trashed anywhere in your Drive.
 * Filters to WorldFerns_ver_*.csv files in the configured folder,
 * then dispatches to GitHub Actions.
 */
function onDriveChange(e) {
  // Only react to newly added files.
  if (e.changeType !== "ADD") return;

  var file;
#### 3c. Paste the script code

1. In the editor you will see a file called `Code.gs` with a default empty
    function. **Select all the existing text and delete it.**
2. Paste the entire block below:

```javascript
// ── Configuration (values are stored in Script Properties, not here) ── //
const GITHUB_TOKEN    = PropertiesService
   .getScriptProperties().getProperty("GITHUB_TOKEN");
const GITHUB_REPO     = "pteridogroup/ppg-import";
const DRIVE_FOLDER_ID = PropertiesService
   .getScriptProperties().getProperty("DRIVE_FOLDER_ID");
// ───────────────────────────────────────────────────────────────────── //

/**
 * Polls the Drive folder for a new WorldFerns file and dispatches to
 * GitHub Actions when one is found. Intended to run on a time-driven
 * trigger every 15 minutes (see step 3f).
 */
function pollDriveAndDispatch() {
   var props  = PropertiesService.getScriptProperties();
   var folder = DriveApp.getFolderById(DRIVE_FOLDER_ID);
   var files  = folder.getFilesByType(MimeType.CSV);

   // Find the newest WorldFerns_ver_*.csv in the folder.
   var newest = null;
   while (files.hasNext()) {
      var f = files.next();
      if (!f.getName().match(/^WorldFerns_ver_.*\.csv$/)) continue;
      if (!newest || f.getLastUpdated() > newest.getLastUpdated()) {
         newest = f;
      }
   }
   if (!newest) {
      Logger.log("No WorldFerns_ver_*.csv found in folder — nothing to do.");
      return;
   }

   // Compare against the last file we dispatched.
   var lastId       = props.getProperty("LAST_DISPATCHED_FILE_ID");
   var lastModified = props.getProperty("LAST_DISPATCHED_MODIFIED");
   var currentModified = newest.getLastUpdated().toISOString();

   if (newest.getId() === lastId && currentModified === lastModified) {
      Logger.log("No new file since last dispatch — skipping.");
      return;
   }

   // Build the GitHub repository_dispatch payload.
   var payload = {
      event_type: "worldferns_updated",
      client_payload: {
         file_id:   newest.getId(),
         file_name: newest.getName(),
         modified:  currentModified
      }
   };

   var resp = UrlFetchApp.fetch(
      "https://api.github.com/repos/" + GITHUB_REPO + "/dispatches",
      {
         method: "post",
         contentType: "application/json",
         headers: {
            Authorization: "Bearer " + GITHUB_TOKEN,
            Accept: "application/vnd.github+json"
         },
         payload: JSON.stringify(payload),
         muteHttpExceptions: true
      }
   );

   if (resp.getResponseCode() >= 300) {
      throw new Error("GitHub dispatch failed: " + resp.getContentText());
   }

   props.setProperty("LAST_DISPATCHED_FILE_ID",  newest.getId());
   props.setProperty("LAST_DISPATCHED_MODIFIED", currentModified);
   Logger.log("Dispatched: " + newest.getName());
}
```

3. Click the **Save** button (disk icon, or Ctrl+S / Cmd+S).

#### 3d. Store secrets in Script Properties

Script Properties are like environment variables — they are stored securely
in your Apps Script project and are never visible in the code.

1. In the Apps Script editor, click the **gear icon (⚙)** in the left sidebar
    to open **Project Settings**.
2. Scroll down to the **Script properties** section and click
    **Add script property**.
3. Add the following two properties one at a time:

    | Property name | Value |
    | ------------- | ----- |
    | `GITHUB_TOKEN` | The GitHub PAT you created in step 2 |
    | `DRIVE_FOLDER_ID` | The folder ID you copied in step 3a |

4. Click **Save script properties**.

> The script will also automatically create `LAST_DISPATCHED_FILE_ID` and
> `LAST_DISPATCHED_MODIFIED` properties on its first successful dispatch —
> you do not need to add those manually.

#### 3e. Authorize the script

Before the trigger can fire, you must grant the script permission to access
Drive and make external HTTP requests.

1. Back in the editor, select `pollDriveAndDispatch` in the function dropdown
    (top toolbar, next to the Run/Debug buttons).
2. Click **Run**.
3. A dialog box will appear: **"Authorization required"**. Click
    **Review permissions**.
4. Choose your Google account.
5. You may see a warning: **"Google hasn't verified this app"**. Click
    **Advanced → Go to (project name) (unsafe)**. This is expected for
    personal scripts — it only means Google has not reviewed the code, not
    that it is dangerous.
6. Click **Allow**.
7. The script will run once. If the folder contains no matching file yet,
    you will see `"No WorldFerns_ver_*.csv found in folder — nothing to
    do."` in the Execution log — that is fine.

#### 3f. Set up the time-driven trigger

1. In the left sidebar, click the **clock icon** (Triggers).
2. Click the **+ Add Trigger** button (bottom right).
3. Fill in the dialog:

    | Setting | Value to select |
    | ------- | --------------- |
    | Choose which function to run | `pollDriveAndDispatch` |
    | Choose which deployment to run | Head |
    | Select event source | **Time-driven** |
    | Select type of time based trigger | **Minutes timer** |
    | Select minute interval | **Every 15 minutes** |
    | Failure notification settings | Notify me daily (recommended) |

4. Click **Save**.

The script will now run every 15 minutes. When a new `WorldFerns_ver_*.csv`
file is uploaded to the Drive folder, the pipeline will start within
15 minutes.

#### 3g. Test the setup end-to-end

1. Upload a `WorldFerns_ver_*.csv` file to the Drive folder.
2. Wait up to 15 minutes for the trigger to fire, or run
    `pollDriveAndDispatch` manually from the editor to trigger immediately.
3. Go to **Apps Script → Executions** (clock icon → Executions tab) to
    confirm the function ran without errors and logged `"Dispatched: ..."`.
4. Go to your GitHub repository → **Actions** tab. You should see a new
    workflow run called **"Import World Ferns from Google Drive"** starting.
5. After it completes, click into the run and download the `wf-dwc-csv`
    artifact to verify the output CSV.

---

## GitHub Actions configuration checklist

Set these in **GitHub → Settings → Secrets and variables → Actions**.

### Repository secrets

- `SHINYAPPS_TOKEN` (required for deploy workflow)
- `SHINYAPPS_SECRET` (required for deploy workflow)
- `GDRIVE_SA_KEY` (required for automated imports from the private Drive
  folder)

### Notes

- `SHINYAPPS_ACCOUNT` and `SHINYAPPS_APP_NAME` are currently hard-coded in
   `.github/workflows/deploy-shinyapps.yml`.
- Apps Script script properties are configured separately (not in GitHub
   Actions): `GITHUB_TOKEN` and `DRIVE_FOLDER_ID`.
