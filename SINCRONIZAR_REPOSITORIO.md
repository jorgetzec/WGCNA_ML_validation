
# Guía Completa: Sincronizar Repositorio Git con Directorio Local (Vinculación Perdida)

## Tabla de Contenidos
1. [Diagnóstico del Problema](#diagnóstico-del-problema)
2. [Configuración y Verificación de Credenciales](#configuración-y-verificación-de-credenciales)
3. [Escenario A: Archivos Idénticos](#escenario-a-archivos-idénticos)
4. [Escenario B: Archivos Diferentes](#escenario-b-archivos-diferentes)
5. [Verificación Final](#verificación-final)
6. [Comandos de Referencia Rápida](#comandos-de-referencia-rápida)

---

## Diagnóstico del Problema

### Paso 1: Verificar si existe un repositorio Git

**Comando:**
```bash
git status
```

**Output esperado si NO hay repositorio:**
```
fatal: not a git repository (or any of the parent directories): .git
```

**Interpretación:**
- Este mensaje indica que Git no encuentra el directorio `.git` en el directorio actual ni en ningún directorio padre.
- El directorio `.git` contiene toda la información del repositorio (historial, ramas, remotos, etc.).
- Si falta este directorio, el repositorio local está desconectado o nunca fue inicializado.

**Comando alternativo (PowerShell):**
```powershell
Test-Path .git
```

**Output esperado:**
- `False` = No existe repositorio Git
- `True` = Existe repositorio Git

### Paso 2: Verificar archivos ocultos

**Comando (PowerShell):**
```powershell
Get-ChildItem -Force -Hidden | Select-Object Name
```

**Interpretación:**
- `-Force` muestra archivos ocultos y del sistema
- `-Hidden` filtra solo archivos ocultos
- Si no aparece `.git` en la lista, confirma que no existe el repositorio

### Paso 3: Identificar el repositorio remoto

**Necesitas conocer:**
- La URL del repositorio remoto (GitHub, GitLab, etc.)
- Ejemplo: `https://github.com/usuario/repositorio.git`

---

## Configuración y Verificación de Credenciales

Cuando trabajas con múltiples cuentas de GitHub (por ejemplo, personal e institucional), es importante saber qué credenciales está usando Git y cómo cambiarlas.

### Verificar Configuración de Git (Usuario y Email)

**Comando:**
```bash
git config --list | grep -E "(user|remote)"
```

**En PowerShell:**
```powershell
git config --list | Select-String -Pattern "(user|remote)"
```

**Output esperado:**
```
user.name=Jorge Tzec
user.email=jorgetzec@gmail.com
remote.origin.url=https://github.com/jorgetzec/WGCNA_ML_validation.git
```

**Interpretación:**
- `user.name` = Nombre que aparecerá en los commits
- `user.email` = Email asociado a los commits
- `remote.origin.url` = URL del repositorio remoto

**Ver solo configuración del repositorio actual:**
```bash
git config user.name
git config user.email
```

**Ver configuración global:**
```bash
git config --global user.name
git config --global user.email
```

### Configurar Credenciales de Git por Repositorio

**Para cuenta principal (Jorge Tzec):**
```bash
git config user.name "Jorge Tzec"
git config user.email "jorgetzec@gmail.com"
```

**Para cuenta secundaria (IPEA):**
```bash
git config user.name "GrupoIPEA"
git config user.email "gestion@grupoipea.com"
```

**Nota:** Estos comandos configuran solo el repositorio actual. Para cambiar globalmente, usa `--global`:
```bash
git config --global user.name "Jorge Tzec"
git config --global user.email "jorgetzec@gmail.com"
```

### Verificar Credenciales HTTPS Guardadas (Windows)

**Comando (PowerShell):**
```powershell
cmdkey /list | Select-String -Pattern "github"
```

**Output esperado:**
```
    Destino: LegacyGeneric:target=git:https://github.com
    Tipo: Genérico
    Usuario: jorgetzec
    Persistencia del equipo local
```

**Interpretación:**
- Muestra qué credenciales de GitHub están guardadas en el Administrador de Credenciales de Windows
- El campo `Usuario` indica qué cuenta está configurada

**Ver todas las credenciales:**
```powershell
cmdkey /list
```

### Eliminar Credenciales HTTPS Guardadas

Si Git está usando credenciales incorrectas (por ejemplo, de otra cuenta), elimínalas:

**Comando (PowerShell):**
```powershell
cmdkey /delete:"LegacyGeneric:target=git:https://github.com"
```

**Después de eliminar:**
- La próxima vez que hagas `git push` o `git pull`, Git pedirá nuevas credenciales
- Ingresa tu usuario de GitHub y un **Personal Access Token (PAT)** como contraseña
- GitHub ya no acepta contraseñas normales para HTTPS

### Configuración SSH para Múltiples Cuentas

Si usas SSH en lugar de HTTPS, puedes configurar múltiples cuentas usando el archivo `~/.ssh/config`:

**Ubicación del archivo (Windows):**
```
C:\Users\[tu_usuario]\.ssh\config
```

**Ejemplo de configuración:**
```
# Cuenta principal (Jorge Tzec)
Host github.com
  HostName github.com
  User git
  IdentityFile ~/.ssh/id_ed25519
  IdentitiesOnly yes

# Cuenta secundaria (IPEA)
Host github-ipea
  HostName github.com
  User git
  IdentityFile ~/.ssh/id_ed25519_ipea
  IdentitiesOnly yes
```

**Usar cuenta principal:**
```bash
git remote set-url origin git@github.com:usuario/repositorio.git
```

**Usar cuenta secundaria:**
```bash
git remote set-url origin git@github-ipea:usuario/repositorio.git
```

**Verificar qué URL está configurada:**
```bash
git remote -v
```

### Verificar Autenticación Actual

**Probar conexión con HTTPS:**
```bash
git ls-remote --heads origin
```

**Probar conexión con SSH:**
```bash
ssh -T git@github.com
```

**Output esperado (SSH):**
```
Hi jorgetzec! You've successfully authenticated, but GitHub does not provide shell access.
```

**Interpretación:**
- Confirma que estás autenticado con la cuenta correcta
- El nombre de usuario en el mensaje indica qué cuenta está activa

### Cambiar entre Cuentas

**Método 1: Cambiar configuración de Git (recomendado para HTTPS)**
```bash
# Configurar para cuenta principal
git config user.name "Jorge Tzec"
git config user.email "jorgetzec@gmail.com"

# Eliminar credenciales guardadas si es necesario
cmdkey /delete:"LegacyGeneric:target=git:https://github.com"

# Hacer push/pull (pedirá nuevas credenciales)
git push origin main
```

**Método 2: Cambiar URL remota (recomendado para SSH)**
```bash
# Para cuenta principal
git remote set-url origin git@github.com:usuario/repositorio.git

# Para cuenta secundaria
git remote set-url origin git@github-ipea:usuario/repositorio.git
```

### Solución de Problemas: Error 403 "Permission denied"

**Síntoma:**
```
remote: Permission to usuario/repositorio.git denied to [cuenta_incorrecta].
fatal: unable to access 'https://github.com/...': The requested URL returned error: 403
```

**Causa:** Git está usando credenciales de otra cuenta.

**Solución:**
1. **Eliminar credenciales guardadas:**
   ```powershell
   cmdkey /delete:"LegacyGeneric:target=git:https://github.com"
   ```

2. **Verificar configuración de Git:**
   ```bash
   git config user.name
   git config user.email
   ```

3. **Configurar correctamente si es necesario:**
   ```bash
   git config user.name "Jorge Tzec"
   git config user.email "jorgetzec@gmail.com"
   ```

4. **Intentar push nuevamente:**
   ```bash
   git push origin main
   ```
   - Git pedirá credenciales nuevas
   - Usa tu usuario de GitHub y un **Personal Access Token** como contraseña

### Resumen: Comandos de Verificación Rápida

```bash
# Ver configuración de usuario/email
git config user.name
git config user.email

# Ver todas las configuraciones
git config --list | grep user

# Ver URL del remoto
git remote -v

# Ver credenciales HTTPS guardadas (Windows)
cmdkey /list | Select-String -Pattern "github"

# Probar autenticación
git ls-remote --heads origin
```

---

## Escenario A: Archivos Idénticos

**Situación:** Los archivos en tu directorio local son exactamente iguales a los del repositorio remoto, pero perdiste la vinculación con Git.

### Proceso Completo

#### Paso 1: Inicializar repositorio Git local

**Comando:**
```bash
git init
```

**Output esperado:**
```
Initialized empty Git repository in [ruta]/[directorio]/.git/
```

**Interpretación:**
- Crea un nuevo directorio `.git` vacío en el directorio actual
- Inicializa un repositorio Git nuevo sin historial
- La rama por defecto será `main` (o `master` en versiones antiguas)
- **Importante:** Esto NO borra tus archivos existentes

**Lógica:**
- Git necesita un directorio `.git` para funcionar
- Este directorio almacena metadatos, historial, configuración, etc.
- Al inicializar, creas la estructura necesaria para que Git pueda rastrear cambios

#### Paso 2: Agregar el repositorio remoto

**Comando:**
```bash
git remote add origin https://github.com/usuario/repositorio.git
```

**Interpretación:**
- `origin` es el nombre estándar para el repositorio remoto principal
- Puedes usar otro nombre, pero `origin` es la convención
- Esto solo configura la URL, NO descarga nada aún

**Verificar que se agregó correctamente:**
```bash
git remote -v
```

**Output esperado:**
```
origin  https://github.com/usuario/repositorio.git (fetch)
origin  https://github.com/usuario/repositorio.git (push)
```

**Interpretación:**
- Muestra todas las URLs remotas configuradas
- `fetch` = URL para descargar cambios
- `push` = URL para subir cambios
- Si ambas son iguales, está correctamente configurado

#### Paso 3: Descargar información del remoto (sin modificar archivos locales)

**Comando:**
```bash
git fetch origin
```

**Output esperado:**
```
From https://github.com/usuario/repositorio
 * [new branch]      main       -> origin/main
 * [new tag]         v1.0.0     -> v1.0.0
```

**Interpretación:**
- `fetch` descarga información del remoto pero NO modifica tus archivos locales
- `origin/main` = referencia remota a la rama `main`
- Los tags también se descargan (ej: `v1.0.0`)
- Tus archivos locales permanecen intactos en este paso

**Lógica:**
- Git separa la descarga de información (`fetch`) de la aplicación de cambios (`merge`/`pull`)
- Esto te permite inspeccionar qué hay en el remoto antes de modificar tu directorio

#### Paso 4: Verificar el estado actual

**Comando:**
```bash
git status
```

**Output esperado:**
```
On branch main

No commits yet

Untracked files:
  (use "git add <file>..." to include in what will be committed)
        archivo1.R
        archivo2.csv
        README.md
        ...

nothing added to commit but untracked files present (use "git add" to track)
```

**Interpretación:**
- `On branch main` = estás en la rama `main` local
- `No commits yet` = no hay historial local aún
- `Untracked files` = archivos que existen pero Git no está rastreando
- Estos archivos son tus archivos locales actuales

**Ver ramas disponibles:**
```bash
git branch -a
```

**Output esperado:**
```
  remotes/origin/HEAD -> origin/main
  remotes/origin/main
```

**Interpretación:**
- `remotes/origin/main` = rama remota descargada con `fetch`
- No hay rama local `main` aún porque no hay commits locales

#### Paso 5: Verificar qué archivos hay en el remoto

**Comando:**
```bash
git ls-tree -r --name-only origin/main
```

**Output esperado:**
```
archivo1.R
archivo2.csv
README.md
...
```

**Interpretación:**
- Lista todos los archivos en la rama `main` del remoto
- `-r` = recursivo (incluye subdirectorios)
- `--name-only` = solo nombres de archivos
- Compara esta lista con tus archivos locales para verificar que son los mismos

**Ver historial del remoto:**
```bash
git log origin/main --oneline -5
```

**Output esperado:**
```
411426e Revise README with citation and version details
d8c2fd6 Adjunt comments in Refactor Random Forest validation script
c592d1b Refine comments and details
1de6d05 Refine comments and messages in analysis script
c73e049 Revise script details in Cre_Salt_GeneOntology.R
```

**Interpretación:**
- Muestra los últimos 5 commits del remoto
- Cada línea = un commit con su hash (código único) y mensaje
- `-5` limita a 5 commits (puedes cambiar el número)

#### Paso 6: Sincronizar con el remoto (archivos idénticos)

**Comando:**
```bash
git reset --hard origin/main
```

**Output esperado:**
```
HEAD is now at 411426e Revise README with citation and version details
```

**Interpretación:**
- `reset --hard` mueve la rama local al mismo commit que `origin/main`
- `--hard` también actualiza el directorio de trabajo para que coincida exactamente
- `HEAD` es el puntero al commit actual
- Ahora tu rama local apunta al mismo commit que el remoto

**Lógica:**
- Como los archivos son idénticos, es seguro hacer un reset hard
- Esto alinea tu historial local con el remoto
- Tus archivos locales se marcan como rastreados por Git

**⚠️ ADVERTENCIA:**
- `git reset --hard` descarta cualquier cambio no guardado
- Solo úsalo cuando estés seguro de que los archivos son idénticos
- Si tienes cambios importantes, haz backup primero

#### Paso 7: Verificar sincronización

**Comando:**
```bash
git status
```

**Output esperado:**
```
On branch main
nothing to commit, working tree clean
```

**Interpretación:**
- `On branch main` = estás en la rama main
- `nothing to commit` = no hay cambios pendientes
- `working tree clean` = tu directorio coincide exactamente con el último commit
- **¡Perfecto!** Estás sincronizado

**Ver historial local:**
```bash
git log --oneline -5
```

**Output esperado:**
```
411426e Revise README with citation and version details
d8c2fd6 Adjunt comments in Refactor Random Forest validation script
c592d1b Refine comments and details
1de6d05 Refine comments and messages in analysis script
c73e049 Revise script details in Cre_Salt_GeneOntology.R
```

**Interpretación:**
- Debe coincidir con el historial del remoto
- Confirma que tu historial local está sincronizado

---

## Escenario B: Archivos Diferentes

**Situación:** Los archivos en tu directorio local son diferentes a los del repositorio remoto, y perdiste la vinculación con Git.

### Opción B1: Preferir archivos remotos (descartar cambios locales)

**⚠️ ADVERTENCIA:** Esto eliminará todos los cambios locales no guardados.

#### Pasos 1-4: Igual que Escenario A
(Seguir pasos 1-4 del Escenario A)

#### Paso 5: Hacer backup de archivos locales (RECOMENDADO)

**Comando (PowerShell):**
```powershell
Copy-Item -Path . -Destination ..\backup_local -Recurse -Exclude .git
```

**Interpretación:**
- Crea una copia de seguridad de tus archivos locales
- `-Exclude .git` evita copiar el directorio Git (no es necesario)
- Ahora puedes recuperar tus archivos si es necesario

#### Paso 6: Eliminar archivos locales que entran en conflicto

**Opción A: Eliminar archivos específicos**
```bash
# Eliminar archivos que están en el remoto
git rm archivo1.R archivo2.csv README.md
```

**Opción B: Eliminar todos los archivos no rastreados**
```bash
# CUIDADO: Esto elimina TODOS los archivos no rastreados
git clean -fd
```

**Interpretación:**
- `git rm` elimina archivos del índice y del sistema de archivos
- `git clean -fd` elimina archivos no rastreados (`-f` = force, `-d` = directorios)
- Esto permite que `git pull` funcione sin conflictos

#### Paso 7: Traer cambios del remoto

**Comando:**
```bash
git pull origin main --allow-unrelated-histories
```

**Output esperado:**
```
From https://github.com/usuario/repositorio
 * branch            main       -> FETCH_HEAD
Merge made by the 'recursive' strategy.
 [archivos nuevos agregados]
```

**Interpretación:**
- `--allow-unrelated-histories` permite fusionar historiales sin ancestro común
- Crea un merge commit que combina ambos historiales
- Los archivos del remoto ahora están en tu directorio local

**Verificar:**
```bash
git status
```

**Output esperado:**
```
On branch main
nothing to commit, working tree clean
```

---

### Opción B2: Preferir archivos locales (guardar cambios locales)

**Situación:** Quieres mantener tus cambios locales y agregarlos al historial.

#### Pasos 1-4: Igual que Escenario A
(Seguir pasos 1-4 del Escenario A)

#### Paso 5: Agregar archivos locales al staging

**Comando:**
```bash
git add .
```

**Interpretación:**
- Agrega todos los archivos al área de staging
- Git ahora rastrea estos archivos
- Preparados para hacer commit

**Verificar qué se agregó:**
```bash
git status
```

**Output esperado:**
```
On branch main

No commits yet

Changes to be committed:
  (use "git rm --cached <file>..." to unstage)
        new file:   archivo1.R
        new file:   archivo2.csv
        new file:   README.md
        ...
```

#### Paso 6: Hacer commit inicial con archivos locales

**Comando:**
```bash
git commit -m "Restaurar archivos locales después de pérdida de vinculación"
```

**Output esperado:**
```
[main (root-commit) abc1234] Restaurar archivos locales después de pérdida de vinculación
 9 files changed, 1234 insertions(+)
 create mode 100644 archivo1.R
 create mode 100644 archivo2.csv
 ...
```

**Interpretación:**
- Crea el primer commit en tu historial local
- `root-commit` indica que es el commit inicial
- `abc1234` es el hash del commit (único)

#### Paso 7: Fusionar con el remoto

**Comando:**
```bash
git pull origin main --allow-unrelated-histories --no-edit
```

**Output esperado:**
```
From https://github.com/usuario/repositorio
 * branch            main       -> FETCH_HEAD
Merge made by the 'recursive' strategy.
 [archivos del remoto agregados]
```

**Interpretación:**
- `--no-edit` usa el mensaje de merge automático (evita abrir editor)
- Crea un merge commit que combina tu historial local con el remoto
- Si hay archivos con el mismo nombre pero contenido diferente, Git intentará fusionarlos automáticamente

#### Paso 8: Resolver conflictos (si los hay)

**Si hay conflictos, verás:**
```
Auto-merging archivo1.R
CONFLICT (add/add): Merge conflict in archivo1.R
Automatic merge failed; fix conflicts and then commit the result.
```

**Proceso de resolución:**

1. **Ver archivos en conflicto:**
```bash
git status
```

**Output:**
```
You have unmerged paths.
  (use "git add <file>..." to mark resolution)

Unmerged paths:
  (use "git add <file>..." to mark resolution)
        both added:      archivo1.R
```

2. **Abrir archivos en conflicto y buscar marcadores:**
```
<<<<<<< HEAD
[tu versión local]
=======
[versión del remoto]
>>>>>>> origin/main
```

3. **Editar manualmente** para decidir qué versión mantener o combinar ambas

4. **Marcar como resuelto:**
```bash
git add archivo1.R
```

5. **Completar el merge:**
```bash
git commit -m "Resolver conflictos de merge"
```

#### Paso 9: Verificar sincronización

**Comando:**
```bash
git status
```

**Output esperado:**
```
On branch main
nothing to commit, working tree clean
```

---

### Opción B3: Comparar y decidir archivo por archivo

**Situación:** Quieres revisar diferencias antes de decidir qué hacer.

#### Pasos 1-4: Igual que Escenario A
(Seguir pasos 1-4 del Escenario A)

#### Paso 5: Comparar archivos específicos

**Ver diferencias de un archivo:**
```bash
git diff --no-index archivo_local.R origin/main:archivo_local.R
```

**O descargar versión remota temporalmente:**
```bash
git show origin/main:archivo.R > archivo_remoto.R
```

Luego comparar manualmente con tu editor.

#### Paso 6: Decidir por archivo

Para cada archivo:
- **Si prefieres remoto:** Eliminar local y hacer pull
- **Si prefieres local:** Agregar al staging y hacer commit
- **Si quieres ambos:** Resolver manualmente después del merge

---

## Verificación Final

### Verificar conexión con remoto

**Comando:**
```bash
git remote show origin
```

**Output esperado:**
```
* remote origin
  Fetch URL: https://github.com/usuario/repositorio.git
  Push  URL: https://github.com/usuario/repositorio.git
  HEAD branch: main
  Remote branch:
    main tracked
  Local ref configured for 'git push':
    main pushes to main (up to date)
```

**Interpretación:**
- `HEAD branch: main` = rama principal del remoto
- `main tracked` = tu rama local rastrea la remota
- `up to date` = estás sincronizado con el remoto

### Verificar que .git existe

**Comando (PowerShell):**
```powershell
Test-Path .git
```

**Output esperado:** `True`

### Verificar estado limpio

**Comando:**
```bash
git status
```

**Output esperado:**
```
On branch main
nothing to commit, working tree clean
```

### Verificar historial

**Comando:**
```bash
git log --oneline --graph -10
```

**Output esperado:**
```
* 411426e Revise README with citation and version details
* d8c2fd6 Adjunt comments in Refactor Random Forest validation script
* c592d1b Refine comments and details
...
```

**Interpretación:**
- `--graph` muestra visualización del historial
- `-10` muestra últimos 10 commits
- Debe coincidir con el remoto (en Escenario A) o incluir tus commits locales (en Escenario B)

---

## Comandos de Referencia Rápida

### Diagnóstico
```bash
git status                          # Ver estado del repositorio
Test-Path .git                      # Verificar si existe .git (PowerShell)
git remote -v                       # Ver remotos configurados
```

### Credenciales
```bash
# Ver configuración de usuario/email
git config user.name                # Ver nombre configurado
git config user.email               # Ver email configurado
git config --list | grep user       # Ver todas las configuraciones de usuario

# Configurar credenciales (repositorio actual)
git config user.name "Jorge Tzec"
git config user.email "jorgetzec@gmail.com"

# Ver credenciales HTTPS guardadas (Windows PowerShell)
cmdkey /list | Select-String -Pattern "github"

# Eliminar credenciales HTTPS guardadas
cmdkey /delete:"LegacyGeneric:target=git:https://github.com"

# Probar autenticación
git ls-remote --heads origin        # Probar conexión con remoto
```

### Inicialización
```bash
git init                            # Inicializar repositorio nuevo
git remote add origin <URL>         # Agregar remoto
git fetch origin                    # Descargar información del remoto
```

### Inspección
```bash
git branch -a                       # Ver todas las ramas (local y remota)
git log origin/main --oneline -5   # Ver historial del remoto
git ls-tree -r --name-only origin/main  # Listar archivos del remoto
```

### Sincronización (Archivos Idénticos)
```bash
git reset --hard origin/main        # Sincronizar con remoto (CUIDADO: descarta cambios)
```

### Sincronización (Archivos Diferentes)
```bash
# Opción 1: Preferir remoto
git clean -fd                       # Eliminar archivos no rastreados
git pull origin main --allow-unrelated-histories

# Opción 2: Preferir local
git add .                           # Agregar archivos locales
git commit -m "Mensaje"             # Hacer commit
git pull origin main --allow-unrelated-histories --no-edit
```

### Verificación
```bash
git status                          # Estado actual
git remote show origin              # Información del remoto
git log --oneline -5                # Historial local
```

---

## Resumen de Decisiones

### ¿Cómo saber qué escenario aplicar?

1. **Ejecuta diagnóstico:**
   ```bash
   git status
   git ls-tree -r --name-only origin/main
   ```

2. **Compara listas de archivos:**
   - Si son idénticos → **Escenario A**
   - Si son diferentes → **Escenario B**

3. **En Escenario B, decide:**
   - ¿Prefieres remoto? → **Opción B1**
   - ¿Prefieres local? → **Opción B2**
   - ¿Quieres revisar? → **Opción B3**

### Tabla de Comandos por Situación

| Situación | Comando Principal | Efecto |
|-----------|-------------------|--------|
| Archivos idénticos | `git reset --hard origin/main` | Sincroniza sin perder nada |
| Preferir remoto | `git clean -fd` + `git pull` | Elimina cambios locales |
| Preferir local | `git add .` + `git commit` + `git pull` | Mantiene cambios locales |
| Revisar diferencias | `git diff` + comparación manual | Control total |

---

## Notas Importantes

1. **Siempre haz backup** antes de usar `git reset --hard` o `git clean -fd`
2. **`--allow-unrelated-histories`** es necesario cuando los historiales no comparten ancestro común
3. **`git fetch`** es seguro: solo descarga información, no modifica archivos
4. **`git pull`** = `git fetch` + `git merge` (combina ambos pasos)
5. El directorio `.git` está oculto por defecto en Windows
6. Los conflictos de merge requieren resolución manual

---

## Solución de Problemas Comunes

### Error: "would be overwritten by merge"
**Causa:** Archivos locales sin rastrear con el mismo nombre que en el remoto
**Solución:** 
```bash
git clean -fd  # Eliminar archivos no rastreados
# O moverlos manualmente antes del pull
```

### Error: "no commit on branch 'main' yet"
**Causa:** Intentas configurar upstream antes de tener commits locales
**Solución:** Primero haz `git reset --hard origin/main` o `git pull`

### Error: "refusing to merge unrelated histories"
**Causa:** Los historiales no comparten ancestro común
**Solución:** Usa `--allow-unrelated-histories` en el comando `git pull`

---

**Última actualización:** Enero 2026
**Autor:** Guía creada para sincronización de repositorios Git
