

def pr(text=""):
    print(text)

def er(text=""):
    print(text)

def warn(msg):
    print(f"⚠️  {msg}")

def fail(msg):
    print(f"❌ ERROR: {msg}")
    raise SystemExit(1)