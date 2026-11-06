/*
 * Presentation-only compound-table search.  The alias dictionary mirrors
 * web/governance_blocked_projection.py and never changes candidate data.
 */
const CompoundSearch = (() => {
  const ZH_COMPOUND_ALIASES = Object.freeze({
    "sinomenine":["青藤碱","清风藤碱"],"magnoflorine":["木兰碱","广玉兰碱"],"berberine":["小檗碱"],"palmatine":["巴马汀"],"jatrorrhizine":["药根碱"],"coptisine":["黄连碱"],"epiberberine":["表小檗碱"],"tetrandrine":["粉防己碱","汉防己甲素"],"fangchinoline":["防己诺林碱"],"dauricine":["蝙蝠葛碱"],"puerarin":["葛根素"],"daidzein":["大豆素"],"daidzin":["大豆苷"],"baicalin":["黄芩苷"],"baicalein":["黄芩素"],"wogonin":["汉黄芩素"],"wogonoside":["汉黄芩苷"],"astragaloside":["黄芪甲苷"],"glycyrrhizic acid":["甘草酸"],"glycyrrhizin":["甘草酸"],"liquiritin":["甘草苷"],"paeoniflorin":["芍药苷"],"albiflorin":["白芍苷"],"aconitine":["乌头碱"],"mesaconitine":["新乌头碱","中乌头碱"],"hypaconitine":["次乌头碱"],"benzoylaconine":["苯甲酰乌头原碱"],"benzoylmesaconine":["苯甲酰新乌头原碱"],"matrine":["苦参碱"],"oxymatrine":["氧化苦参碱"],"sophocarpine":["槐果碱"],"ginsenoside":["人参皂苷"],"notoginsenoside":["三七皂苷"],"salidroside":["红景天苷"],"tetramethylpyrazine":["川芎嗪"],"ligustrazine":["川芎嗪"],"ferulic acid":["阿魏酸"],"tetrahydropalmatine":["延胡索乙素"],"corydaline":["延胡索甲素","紫堇碱"],"cryptotanshinone":["隐丹参酮"],"tanshinone":["丹参酮"],"danshensu":["丹参素"],"salvianolic acid":["丹酚酸"],"protocatechuic acid":["原儿茶酸"],"hesperidin":["橙皮苷"],"hesperetin":["橙皮素"],"naringin":["柚皮苷"],"neohesperidin":["新橙皮苷"],"synephrine":["辛弗林"],"ephedrine":["麻黄碱"],"pseudoephedrine":["伪麻黄碱"],"quercetin":["槲皮素"],"kaempferol":["山奈酚"],"rutin":["芦丁"],"hyperoside":["金丝桃苷"],"isoquercitrin":["异槲皮苷"],"vitexin":["牡荆素"],"emodin":["大黄素"],"aloe-emodin":["芦荟大黄素"],"rhein":["大黄酸"],"chrysophanol":["大黄酚"],"physcion":["大黄素甲醚"],"ligustilide":["藁本内酯"],"catalpol":["梓醇"],"acteoside":["毛蕊花糖苷"],"verbascoside":["毛蕊花糖苷"],"echinacoside":["松果菊苷"],"forsythiaside":["连翘酯苷"],"phillyrin":["连翘苷"],"arctiin":["牛蒡苷"],"arctigenin":["牛蒡苷元"],"chlorogenic acid":["绿原酸"],"caffeic acid":["咖啡酸"],"gallic acid":["没食子酸"],"ellagic acid":["鞣花酸"],"bergenin":["岩白菜素"],"schisandrin":["五味子素"],"deoxyschisandrin":["五味子甲素"],"curcumin":["姜黄素"],"gingerol":["姜辣素"],"aucubin":["桃叶珊瑚苷"],"geniposide":["栀子苷"],"genipin":["栀子苷元"],"ursolic acid":["熊果酸"],"oleanolic acid":["齐墩果酸"],"betulinic acid":["白桦脂酸"],"sitosterol":["谷甾醇"],"stigmasterol":["豆甾醇"]
  });

  const text = value => String(value || '').trim().toLocaleLowerCase();

  function aliasesFor(names) {
    const aliases = [];
    for (const name of names || []) {
      const normalized = text(name);
      if (!normalized) continue;
      const tokens = new Set(normalized.split(/[^a-z0-9]+/).filter(Boolean));
      tokens.add(normalized);
      for (const key of Object.keys(ZH_COMPOUND_ALIASES).sort((a, b) => b.length - a.length)) {
        if (!tokens.has(key)) continue;
        aliases.push(...ZH_COMPOUND_ALIASES[key]);
        tokens.delete(key);
      }
    }
    return aliases;
  }

  function matches({ names = [], terms = [], mzs = [] } = {}, query) {
    const needle = text(query);
    if (!needle) return true;
    const haystack = [...names, ...terms, ...aliasesFor(names)].map(text).filter(Boolean);
    if (haystack.some(value => value.includes(needle))) return true;
    const numeric = needle.match(/^(?:m\/z\s*)?(\d+(?:\.\d+)?)$/i);
    if (!numeric) return false;
    const target = Number(numeric[1]);
    return mzs.some(value => Number.isFinite(Number(value)) && Math.abs(Number(value) - target) <= 0.5);
  }

  return { ZH_COMPOUND_ALIASES, matches };
})();
