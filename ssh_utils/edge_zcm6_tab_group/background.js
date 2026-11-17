// pcopen (zcm6 -> ssh pc) opens https://zcm6.invalid/open?u=<encoded url> instead of the real URL.
// A marker URL, not a suffix on the real one, so no site's own routing sees anything extra.
// This worker catches it, sends the tab to the real URL and moves it into the window's "zcm6" group.
const MARKER = "https://zcm6.invalid/open?";
const GROUP_TITLE = "zcm6";

async function groupTab(tab) {
  const [group] = await chrome.tabGroups.query({ windowId: tab.windowId, title: GROUP_TITLE });
  if (group) {
    await chrome.tabs.group({ groupId: group.id, tabIds: tab.id });
  } else {
    const groupId = await chrome.tabs.group({ tabIds: tab.id, createProperties: { windowId: tab.windowId } });
    await chrome.tabGroups.update(groupId, { title: GROUP_TITLE, color: "cyan" });
  }
}

chrome.tabs.onUpdated.addListener(async (tabId, changeInfo, tab) => {
  const url = changeInfo.url || "";
  if (!url.startsWith(MARKER)) return;
  const target = new URLSearchParams(url.slice(MARKER.length)).get("u");
  // Only forward real web pages; anything else stays on the dead marker page.
  if (!target || !/^https?:\/\//i.test(target)) return;
  await chrome.tabs.update(tabId, { url: target });
  await groupTab(tab);
});
