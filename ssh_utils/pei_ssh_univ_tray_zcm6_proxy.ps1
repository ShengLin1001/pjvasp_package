# zcm6-proxy 托盘：用隐藏的 mintty 窗口跑 pei_ssh_univ_keep_zcm6_proxy，只在右下角通知区（折叠栏）留一个图标。
# 左键单击图标显示/隐藏日志窗口；右键菜单"退出"关闭窗口即停止隧道。窗口被关掉（含手动点 X）时图标自动消失。
# 已有 zcm6-proxy 窗口（如旧版最小化启动的）就直接接管并藏起来，不再开第二个。
# 开机自启：shell:startup\zcm6-proxy.lnk → powershell -NoProfile -WindowStyle Hidden -ExecutionPolicy Bypass -File <本文件>
# 必须存成 UTF-8 with BOM：Windows PowerShell 5.1 把无 BOM 的文件按 GBK 读，中文菜单会乱码

Add-Type -AssemblyName System.Windows.Forms, System.Drawing
Add-Type -Namespace Win32 -Name User32 -MemberDefinition @'
[DllImport("user32.dll", CharSet = CharSet.Unicode)] public static extern IntPtr FindWindow(string cls, string title);
[DllImport("user32.dll")] public static extern bool IsWindow(IntPtr hwnd);
[DllImport("user32.dll")] public static extern bool IsWindowVisible(IntPtr hwnd);
[DllImport("user32.dll")] public static extern bool IsIconic(IntPtr hwnd);
[DllImport("user32.dll")] public static extern bool ShowWindow(IntPtr hwnd, int cmd);
[DllImport("user32.dll")] public static extern bool SetForegroundWindow(IntPtr hwnd);
[DllImport("user32.dll")] public static extern bool PostMessage(IntPtr hwnd, uint msg, IntPtr w, IntPtr l);
'@

$TITLE = 'zcm6-proxy'
$MINTTY = 'D:\Program Files\Git\usr\bin\mintty.exe'
$ICON = Join-Path $PSScriptRoot 'zcm6-proxy.ico'
# ConfirmExit=no：窗口隐藏时不弹"仍有进程在运行"确认框；关窗口时 mintty 给 bash/ssh 发 SIGHUP，隧道随之断开
$MINTTY_ARGS = "--hold error --window hide -o ConfirmExit=no --title $TITLE -i `"$ICON`" /usr/bin/bash -lc pei_ssh_univ_keep_zcm6_proxy"

function get_window {
    # mintty 刚启动时窗口还没建好，最多等 10s
    for ($i = 0; $i -lt 50; $i++) {
        $hwnd = [Win32.User32]::FindWindow('mintty', $TITLE)
        if ($hwnd -ne [IntPtr]::Zero) { return $hwnd }
        Start-Sleep -Milliseconds 200
    }
    return [IntPtr]::Zero
}

function switch_window {
    if ([Win32.User32]::IsWindowVisible($hwnd) -and -not [Win32.User32]::IsIconic($hwnd)) {
        [void][Win32.User32]::ShowWindow($hwnd, 0)    # SW_HIDE：连任务栏按钮一起藏
    } else {
        [void][Win32.User32]::ShowWindow($hwnd, 9)    # SW_RESTORE：隐藏或最小化都能恢复
        [void][Win32.User32]::SetForegroundWindow($hwnd)
    }
}

### check: 托盘单实例，重复启动直接退出，避免两个图标
$mutex = New-Object System.Threading.Mutex($false, 'Local\pei_ssh_univ_tray_zcm6_proxy')
if (-not $mutex.WaitOne(0)) { exit 0 }
### to here

### prepare: 接管已有窗口，没有才新开
$hwnd = [Win32.User32]::FindWindow('mintty', $TITLE)
if ($hwnd -eq [IntPtr]::Zero) {
    Start-Process $MINTTY -ArgumentList $MINTTY_ARGS
    $hwnd = get_window
}
if ($hwnd -eq [IntPtr]::Zero) {
    [void][System.Windows.Forms.MessageBox]::Show("❌ 10s 内没等到 $TITLE 窗口，检查 $MINTTY 是否存在", $TITLE)
    exit 1
}
[void][Win32.User32]::ShowWindow($hwnd, 0)
### to here

### main: 托盘图标 + 菜单，窗口没了就退出
$tray = New-Object System.Windows.Forms.NotifyIcon
$tray.Icon = New-Object System.Drawing.Icon $ICON
$tray.Text = "$TITLE 反向隧道（单击显示/隐藏窗口）"
$menu = New-Object System.Windows.Forms.ContextMenuStrip
[void]$menu.Items.Add('显示/隐藏窗口', $null, { switch_window })
[void]$menu.Items.Add('退出（停止隧道）', $null, { [void][Win32.User32]::PostMessage($hwnd, 0x10, [IntPtr]::Zero, [IntPtr]::Zero) })  # WM_CLOSE
$tray.ContextMenuStrip = $menu
$tray.add_MouseClick({ if ($_.Button -eq 'Left') { switch_window } })
$tray.Visible = $true

$timer = New-Object System.Windows.Forms.Timer
$timer.Interval = 2000
$timer.add_Tick({ if (-not [Win32.User32]::IsWindow($hwnd)) { [System.Windows.Forms.Application]::Exit() } })
$timer.Start()
[System.Windows.Forms.Application]::Run()
$tray.Dispose()
### to here
