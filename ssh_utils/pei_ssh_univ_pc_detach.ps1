# Run a command on this PC so it outlives the ssh session that started it (zcm6 -> `ssh pc`).
# Why: Win32-OpenSSH puts every session in a job object that kills all children on disconnect, so
# `ssh pc start ...` launches OVITO/Edge and they die as soon as ssh returns. CREATE_BREAKAWAY_FROM_JOB
# moves the new process out of that job; the sshd job allows breakaway (claude --bg survives too).
# Usage (from zcm6):
#   ssh pc 'powershell -NoProfile -File F:\...\pei_ssh_univ_pc_detach.ps1 D:\software\Ovito\ovito-zcm6.cmd sftp://zcm6-ovito/<abs path>'
#   ssh pc 'powershell -NoProfile -File F:\...\pei_ssh_univ_pc_detach.ps1 cmd /c start "" https://example.org'
# The command runs via `cmd /c` with no console window; GUI programs it starts show on the desktop as usual.
# Messages stay ASCII: Windows PowerShell 5.1 reads a BOM-less .ps1 as ANSI and mangles emoji.
if ($args.Count -eq 0) { Write-Host "ERROR: no command given"; exit 1 }

Add-Type @'
using System; using System.Runtime.InteropServices;
public static class PcDetach {
    [StructLayout(LayoutKind.Sequential, CharSet = CharSet.Unicode)]
    struct STARTUPINFO { public int cb; public string r, d, t; public int x, y, w, h, cx, cy, fill, flags; public short show, r2; public IntPtr r3, i, o, e; }
    [StructLayout(LayoutKind.Sequential)]
    struct PROCESS_INFORMATION { public IntPtr hProcess, hThread; public int pid, tid; }
    [DllImport("kernel32.dll", SetLastError = true, CharSet = CharSet.Unicode)]
    static extern bool CreateProcess(string app, string cmd, IntPtr pa, IntPtr ta, bool inherit, uint flags, IntPtr env, string cwd, ref STARTUPINFO si, out PROCESS_INFORMATION pi);
    [DllImport("kernel32.dll")] static extern bool CloseHandle(IntPtr h);
    const uint CREATE_BREAKAWAY_FROM_JOB = 0x01000000, CREATE_NO_WINDOW = 0x08000000;
    public static int Start(string cmdline) {
        var si = new STARTUPINFO(); si.cb = Marshal.SizeOf(si);
        PROCESS_INFORMATION pi;
        if (!CreateProcess(null, cmdline, IntPtr.Zero, IntPtr.Zero, false, CREATE_BREAKAWAY_FROM_JOB | CREATE_NO_WINDOW, IntPtr.Zero, null, ref si, out pi))
            throw new System.ComponentModel.Win32Exception(Marshal.GetLastWin32Error());
        CloseHandle(pi.hThread); CloseHandle(pi.hProcess);
        return pi.pid;
    }
}
'@

# Args arrive already split by PowerShell; requote any with spaces or cmd metacharacters (URLs carry &) before handing them to cmd.
$cmdline = ($args | ForEach-Object { if ("$_" -match '[\s&|<>^]') { '"' + $_ + '"' } else { "$_" } }) -join ' '
$cmdline = "cmd.exe /c " + $cmdline
try {
    $procId = [PcDetach]::Start($cmdline)
    Write-Host "OK: started outside the ssh job (PID $procId): $cmdline"
} catch {
    Write-Host "ERROR: CreateProcess failed: $($_.Exception.Message)"; exit 1
}
